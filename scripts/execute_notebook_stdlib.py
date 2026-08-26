"""Execute this repository's simple Python notebook without Jupyter dependencies."""

from __future__ import annotations

import ast
import contextlib
import io
import json
import os
import traceback
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
NOTEBOOK = ROOT / "notebooks" / "final_signed_six_reproduction.ipynb"


def source_text(cell: dict[str, object]) -> str:
    source = cell.get("source", "")
    return "".join(source) if isinstance(source, list) else str(source)


def main() -> int:
    notebook = json.loads(NOTEBOOK.read_text(encoding="utf-8"))
    if notebook.get("nbformat") != 4 or not isinstance(notebook.get("cells"), list):
        raise ValueError("Invalid notebook structure")

    namespace: dict[str, object] = {"__name__": "notebook_execution"}
    execution_count = 0
    original_cwd = Path.cwd()
    os.chdir(ROOT)
    try:
        for cell in notebook["cells"]:
            if cell.get("cell_type") != "code":
                continue
            execution_count += 1
            source = source_text(cell)
            output_buffer = io.StringIO()
            outputs: list[dict[str, object]] = []
            try:
                tree = ast.parse(source, filename=str(NOTEBOOK))
                final_expression = None
                if tree.body and isinstance(tree.body[-1], ast.Expr):
                    final_expression = ast.Expression(tree.body.pop().value)
                with contextlib.redirect_stdout(output_buffer), contextlib.redirect_stderr(output_buffer):
                    if tree.body:
                        exec(compile(tree, str(NOTEBOOK), "exec"), namespace)
                    value = (
                        eval(compile(final_expression, str(NOTEBOOK), "eval"), namespace)
                        if final_expression is not None
                        else None
                    )
                stream = output_buffer.getvalue()
                if stream:
                    outputs.append({"name": "stdout", "output_type": "stream", "text": stream})
                if value is not None:
                    outputs.append(
                        {
                            "data": {"text/plain": repr(value)},
                            "execution_count": execution_count,
                            "metadata": {},
                            "output_type": "execute_result",
                        }
                    )
            except Exception as exc:  # noqa: BLE001 - notebook must preserve the exact failure
                outputs.append(
                    {
                        "ename": type(exc).__name__,
                        "evalue": str(exc),
                        "output_type": "error",
                        "traceback": traceback.format_exc().splitlines(),
                    }
                )
                cell["outputs"] = outputs
                cell["execution_count"] = execution_count
                NOTEBOOK.write_text(
                    json.dumps(notebook, ensure_ascii=False, indent=1) + "\n",
                    encoding="utf-8",
                    newline="\n",
                )
                raise
            cell["outputs"] = outputs
            cell["execution_count"] = execution_count
    finally:
        os.chdir(original_cwd)

    encoded = json.dumps(notebook, ensure_ascii=False, indent=1) + "\n"
    private_markers = ["C:" + "\\\\" + "Users" + "\\\\", "Zyw" + "01"]
    if any(marker in encoded for marker in private_markers):
        raise ValueError("Notebook output contains a private absolute path")
    NOTEBOOK.write_text(encoded, encoding="utf-8", newline="\n")
    print(f"NOTEBOOK EXECUTION: PASS | code_cells={execution_count}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
