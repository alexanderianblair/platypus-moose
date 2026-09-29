#!/usr/bin/env python3
# * This file is part of the MOOSE framework
# * https://mooseframework.inl.gov
# *
# * All rights reserved, see COPYRIGHT for full restrictions
# * https://github.com/idaholab/moose/blob/master/COPYRIGHT
# *
# * Licensed under LGPL 2.1, please see LICENSE for details
# * https://www.gnu.org/licenses/lgpl-2.1.html
"""Structural check of JSON-encoded VTK cell-grid (.dg) output from MFEMCellGridDataCollection.

Uses only the Python standard library. Checks the composite index, that every array has
consistent sizes and finite values, that connectivity indices are valid, that every attribute is
defined on every cell type and refers to existing arrays, and optional expectations such as the
total number of cells of a type or the function space used for a variable.

Usage:
  check_cellgrid.py INDEX.dg [--cells TYPE=N] [--space VAR=SPACE:ORDER] ...
"""
import argparse, json, math, os, sys


def fail(msg):
    print("FAIL: " + msg)
    sys.exit(1)


def load_json(path):
    with open(path, "rb") as f:
        head = f.read(1)
    if head != b"{":
        fail(f"{path} is not JSON (use encoding = JSON for this check)")
    with open(path) as f:
        return json.load(f)


def check_leaf(path, doc):
    if doc.get("data-type") != "cell-grid":
        fail(f"{path}: data-type is {doc.get('data-type')!r}")
    arrays = {}
    for group, lst in doc["arrays"].items():
        for a in lst:
            n, comps, tuples = len(a["data"]), a["components"], a["tuples"]
            if n != comps * tuples:
                fail(f"{path}: array {group}/{a['name']} has {n} values, expected {comps}*{tuples}")
            if a["type"] == "double" and not all(math.isfinite(v) for v in a["data"]):
                fail(f"{path}: array {group}/{a['name']} has non-finite values")
            arrays[(group, a["name"])] = a
    npoints = arrays[("points", "coords")]["tuples"] if ("points", "coords") in arrays else 0
    types = {}
    for ct in doc["cell-types"]:
        conn = arrays.get(tuple(ct["cell-spec"]["connectivity"]))
        if conn is None:
            fail(f"{path}: missing connectivity for {ct['type']}")
        if conn["data"] and (min(conn["data"]) < 0 or max(conn["data"]) >= npoints):
            fail(f"{path}: connectivity of {ct['type']} out of range [0, {npoints})")
        types[ct["type"]] = conn["tuples"]
    spaces = {}
    shape_count = 0
    for att in doc["attributes"]:
        shape_count += bool(att.get("shape"))
        for ct in types:
            info = att["cell-info"].get(ct)
            if info is None:
                fail(f"{path}: attribute {att['name']} is not defined on {ct}")
            for role, ref in info.get("arrays", {}).items():
                arr = arrays.get(tuple(ref))
                if arr is None:
                    fail(f"{path}: attribute {att['name']} role {role} refers to missing {ref}")
                if role == "values" and "dof-sharing" not in info and arr["tuples"] != types[ct]:
                    fail(f"{path}: attribute {att['name']} has {arr['tuples']} tuples on {ct}, "
                         f"expected {types[ct]}")
            spaces.setdefault(att["name"], set()).add(
                f"{info['function-space']}:{info['order']}")
    if shape_count != 1:
        fail(f"{path}: expected exactly one shape attribute, found {shape_count}")
    return types, spaces


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("index")
    ap.add_argument("--cells", action="append", default=[], help="TYPE=N total cells of TYPE")
    ap.add_argument("--space", action="append", default=[], help="VAR=SPACE:ORDER")
    args = ap.parse_args()

    top = load_json(args.index)
    root = os.path.dirname(os.path.abspath(args.index))
    if top.get("data-type") == "composite":
        files = [os.path.join(root, f) for f in top["group"]["files"]]
        if not files:
            fail("composite index lists no pieces")
    else:
        files = [args.index]

    total, spaces = {}, {}
    for f in files:
        if not os.path.exists(f):
            fail(f"piece {f} does not exist")
        types, sp = check_leaf(f, load_json(f))
        for k, v in types.items():
            total[k] = total.get(k, 0) + v
        for k, v in sp.items():
            spaces.setdefault(k, set()).update(v)

    for spec in args.cells:
        ct, n = spec.split("=")
        if total.get(ct, 0) != int(n):
            fail(f"expected {n} cells of type {ct}, found {total.get(ct, 0)}")
    for spec in args.space:
        var, want = spec.split("=")
        if spaces.get(var) != {want}:
            fail(f"variable {var}: expected function space {want}, found {spaces.get(var)}")
    print(f"OK: {len(files)} piece(s), cells {total}, attributes {sorted(spaces)}")


if __name__ == "__main__":
    main()
