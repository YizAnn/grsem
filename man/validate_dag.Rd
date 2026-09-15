\name{validate_dag}
\alias{validate_dag}
\title{Diagnose directed cycles in a structural regression graph}
\usage{validate_dag(edges, nodes = character())}
\arguments{
\item{edges}{Data frame with character from and to columns. Node names must
be nonempty and nonmissing. Duplicate edges are treated as a single edge.}
\item{nodes}{Optional extra nodes, including isolated vertices.}
}
\description{Checks acyclicity using topological sorting and extracts a concrete
directed cycle with iterative depth-first search when necessary.}
\value{List containing is_dag, topological_order, cycle and nodes.
Non-DAGs have an empty topological_order and a closed cycle with the start node
repeated at the end. DAGs have an empty cycle. The empty graph is a DAG.}
\details{This is legality diagnosis, not NOTEARS, structure learning, causal
discovery, or proof of statistical identification.}
