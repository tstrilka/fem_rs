use anyhow::{Context, Result, anyhow, bail};
use std::collections::{HashMap, HashSet};
use std::path::Path;

use crate::mesh_2d::Mesh2D;

/// gmsh element type 2 = 3-node triangle.
const ELEM_TRIANGLE: u32 = 2;

/// Read a gmsh `.msh` v4.1 ASCII file and return a 2D triangular mesh.
///
/// Only triangle elements (type 2) are loaded — lines, points, and
/// higher-order elements are ignored.
/// Z coordinates are dropped (we assume planar 2D meshes).
/// Boundary nodes are computed topologically: any edge appearing in
/// exactly one triangle is a boundary edge.
pub fn parse_msh<P: AsRef<Path>>(path: P) -> Result<Mesh2D> {
    let path = path.as_ref();
    let content = std::fs::read_to_string(path)
        .with_context(|| format!("failed to read {path:?}"))?;
    parse_msh_str(&content)
}

pub fn parse_msh_str(content: &str) -> Result<Mesh2D> {
    let mut cur = Cursor::new(content);

    let mut nodes_by_tag: HashMap<usize, [f64; 2]> = HashMap::new();
    let mut triangle_tags: Vec<[usize; 3]> = Vec::new();
    let mut format_seen = false;

    while let Some(line) = cur.peek() {
        let trimmed = line.trim();
        match trimmed {
            "$MeshFormat" => {
                cur.advance();
                parse_mesh_format(&mut cur)?;
                format_seen = true;
            }
            "$Nodes" => {
                cur.advance();
                parse_nodes(&mut cur, &mut nodes_by_tag)?;
            }
            "$Elements" => {
                cur.advance();
                parse_elements(&mut cur, &mut triangle_tags)?;
            }
            s if s.starts_with('$') && !s.starts_with("$End") => {
                cur.advance();
                let end_marker = format!("$End{}", &s[1..]);
                cur.skip_to(&end_marker)?;
            }
            _ => {
                cur.advance();
            }
        }
    }

    if !format_seen {
        bail!("no $MeshFormat section found");
    }
    if nodes_by_tag.is_empty() {
        bail!("no nodes parsed");
    }
    if triangle_tags.is_empty() {
        bail!("no triangle elements (type 2) parsed");
    }

    let mut sorted_tags: Vec<usize> = nodes_by_tag.keys().copied().collect();
    sorted_tags.sort_unstable();
    let tag_to_idx: HashMap<usize, usize> = sorted_tags
        .iter()
        .enumerate()
        .map(|(i, &t)| (t, i))
        .collect();
    let nodes: Vec<[f64; 2]> = sorted_tags.iter().map(|t| nodes_by_tag[t]).collect();

    let triangles: Vec<[usize; 3]> = triangle_tags
        .iter()
        .map(|t| {
            let to_idx = |tag: usize| -> Result<usize> {
                tag_to_idx
                    .get(&tag)
                    .copied()
                    .ok_or_else(|| anyhow!("triangle references unknown node tag {tag}"))
            };
            Ok([to_idx(t[0])?, to_idx(t[1])?, to_idx(t[2])?])
        })
        .collect::<Result<Vec<_>>>()?;

    let boundary_nodes = compute_boundary_nodes(&triangles);

    Ok(Mesh2D {
        nodes,
        triangles,
        boundary_nodes,
    })
}

fn parse_mesh_format(cur: &mut Cursor) -> Result<()> {
    let header = cur.next()?;
    let parts: Vec<&str> = header.split_whitespace().collect();
    if parts.len() < 3 {
        bail!("malformed $MeshFormat header: {header}");
    }
    let version = parts[0];
    if !version.starts_with("4.") {
        bail!("only gmsh format 4.x supported, got version {version}");
    }
    let file_type: u32 = parts[1]
        .parse()
        .with_context(|| format!("file-type field: {}", parts[1]))?;
    if file_type != 0 {
        bail!("only ASCII format supported (file-type=0), got file-type={file_type}");
    }
    cur.skip_to("$EndMeshFormat")
}

fn parse_nodes(cur: &mut Cursor, out: &mut HashMap<usize, [f64; 2]>) -> Result<()> {
    let header = cur.next()?;
    let h: Vec<&str> = header.split_whitespace().collect();
    if h.len() < 4 {
        bail!("malformed $Nodes header: {header}");
    }
    let num_blocks: usize = h[0]
        .parse()
        .with_context(|| format!("$Nodes block count: {}", h[0]))?;

    for _ in 0..num_blocks {
        let block_hdr = cur.next()?;
        let b: Vec<&str> = block_hdr.split_whitespace().collect();
        if b.len() < 4 {
            bail!("malformed node block header: {block_hdr}");
        }
        let num_in_block: usize = b[3]
            .parse()
            .with_context(|| format!("node block count: {}", b[3]))?;

        let mut tags = Vec::with_capacity(num_in_block);
        for _ in 0..num_in_block {
            let l = cur.next()?;
            let tag: usize = l
                .trim()
                .parse()
                .with_context(|| format!("node tag: {l}"))?;
            tags.push(tag);
        }
        for tag in tags {
            let l = cur.next()?;
            let c: Vec<&str> = l.split_whitespace().collect();
            if c.len() < 3 {
                bail!("malformed node coordinates: {l}");
            }
            let x: f64 = c[0].parse().with_context(|| format!("node x: {}", c[0]))?;
            let y: f64 = c[1].parse().with_context(|| format!("node y: {}", c[1]))?;
            out.insert(tag, [x, y]);
        }
    }
    cur.skip_to("$EndNodes")
}

fn parse_elements(cur: &mut Cursor, out: &mut Vec<[usize; 3]>) -> Result<()> {
    let header = cur.next()?;
    let h: Vec<&str> = header.split_whitespace().collect();
    if h.len() < 4 {
        bail!("malformed $Elements header: {header}");
    }
    let num_blocks: usize = h[0]
        .parse()
        .with_context(|| format!("$Elements block count: {}", h[0]))?;

    for _ in 0..num_blocks {
        let block_hdr = cur.next()?;
        let b: Vec<&str> = block_hdr.split_whitespace().collect();
        if b.len() < 4 {
            bail!("malformed element block header: {block_hdr}");
        }
        let elem_type: u32 = b[2]
            .parse()
            .with_context(|| format!("element type: {}", b[2]))?;
        let num_in_block: usize = b[3]
            .parse()
            .with_context(|| format!("element block count: {}", b[3]))?;

        for _ in 0..num_in_block {
            let l = cur.next()?;
            if elem_type == ELEM_TRIANGLE {
                let parts: Vec<&str> = l.split_whitespace().collect();
                if parts.len() < 4 {
                    bail!("malformed triangle: {l}");
                }
                let n0: usize = parts[1].parse()?;
                let n1: usize = parts[2].parse()?;
                let n2: usize = parts[3].parse()?;
                out.push([n0, n1, n2]);
            }
        }
    }
    cur.skip_to("$EndElements")
}

fn compute_boundary_nodes(triangles: &[[usize; 3]]) -> Vec<usize> {
    let mut edge_count: HashMap<(usize, usize), u32> = HashMap::new();
    for tri in triangles {
        for k in 0..3 {
            let a = tri[k];
            let b = tri[(k + 1) % 3];
            let key = if a < b { (a, b) } else { (b, a) };
            *edge_count.entry(key).or_insert(0) += 1;
        }
    }
    let mut nodes_set: HashSet<usize> = HashSet::new();
    for ((a, b), count) in edge_count {
        if count == 1 {
            nodes_set.insert(a);
            nodes_set.insert(b);
        }
    }
    let mut nodes: Vec<usize> = nodes_set.into_iter().collect();
    nodes.sort_unstable();
    nodes
}

struct Cursor<'a> {
    lines: Vec<&'a str>,
    i: usize,
}

impl<'a> Cursor<'a> {
    fn new(s: &'a str) -> Self {
        Self {
            lines: s.lines().collect(),
            i: 0,
        }
    }

    fn peek(&self) -> Option<&'a str> {
        self.lines.get(self.i).copied()
    }

    fn advance(&mut self) {
        self.i += 1;
    }

    fn next(&mut self) -> Result<&'a str> {
        let l = self
            .lines
            .get(self.i)
            .copied()
            .ok_or_else(|| anyhow!("unexpected end of file"))?;
        self.i += 1;
        Ok(l)
    }

    fn skip_to(&mut self, marker: &str) -> Result<()> {
        while let Ok(l) = self.next() {
            if l.trim() == marker {
                return Ok(());
            }
        }
        bail!("missing closing marker {marker}");
    }
}
