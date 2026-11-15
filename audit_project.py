
from pathlib import Path
import ast
import csv
import subprocess
import zipfile
from datetime import datetime

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "project_audit"
OUT.mkdir(exist_ok=True)

SKIP = {
    ".git", ".venv", "venv", "__pycache__",
    "node_modules", ".idea", ".vscode",
    "phyloresist", "site-packages"
}

def run_git(*args):
    try:
        return subprocess.run(
            ["git", *args],
            cwd=ROOT,
            capture_output=True,
            text=True,
            timeout=15
        ).stdout.strip()
    except Exception:
        return "Git information unavailable."

def scan_files():
    files = []
    for p in ROOT.rglob("*"):
        if any(part in SKIP for part in p.relative_to(ROOT).parts):
            continue
        if p.is_file() and OUT not in p.parents:
            files.append(p)
    return sorted(files)

files = scan_files()
lines = [
    "# PhyloResist AI - Project Audit",
    f"Generated: {datetime.now().isoformat()}",
    f"Root: {ROOT}",
    f"Total files scanned: {len(files)}",
    "",
    "## 1. Repository structure",
]

for p in files:
    rel = p.relative_to(ROOT)
    lines.append(f"- {rel} ({p.stat().st_size:,} bytes)")

lines += ["", "## 2. Python code architecture"]

for p in files:
    if p.suffix != ".py":
        continue
    try:
        source = p.read_text(encoding="utf-8")
        tree = ast.parse(source)
        rel = p.relative_to(ROOT)
        lines.append(f"\n### {rel}")
        for node in tree.body:
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                lines.append(f"- Function: {node.name} (line {node.lineno})")
            elif isinstance(node, ast.ClassDef):
                lines.append(f"- Class: {node.name} (line {node.lineno})")
            elif isinstance(node, (ast.Import, ast.ImportFrom)):
                lines.append(
                    f"- Import: {ast.get_source_segment(source, node)}"
                )
    except Exception as e:
        lines.append(f"- Could not parse {p.name}: {e}")

lines += ["", "## 3. Dataset inventory (metadata only)"]

for p in files:
    if p.suffix.lower() == ".csv":
        try:
            with p.open(encoding="utf-8-sig", newline="") as f:
                reader = csv.reader(f)
                header = next(reader, [])
            lines.append(
                f"- {p.relative_to(ROOT)}: columns={header}"
            )
        except Exception as e:
            lines.append(f"- {p.name}: inspection error: {e}")

lines += ["", "## 4. Project documentation"]

total_doc_chars = 0
for p in files:
    if p.name.lower() in {"readme.md", "readme.txt"}:
        try:
            content = p.read_text(encoding="utf-8")
            lines.append(f"\n### {p.relative_to(ROOT)}")
            lines.append(content[:12000])
        except Exception:
            pass

for p in files:
    if p.suffix.lower() != ".md":
        continue
    if "results" not in p.parts and "reports" not in p.parts:
        continue
    try:
        content = p.read_text(encoding="utf-8")
        remaining = max(0, 60000 - total_doc_chars)
        if remaining == 0:
            break
        content = content[:min(15000, remaining)]
        lines.append(f"\n### Report: {p.relative_to(ROOT)}")
        lines.append(content)
        total_doc_chars += len(content)
    except Exception:
        pass

lines += [
    "",
    "## 5. Git status",
    run_git("status", "--short"),
    "",
    "## 6. Recent commits",
    run_git("log", "-10", "--oneline"),
    "",
    "## 7. Recent change summary",
    run_git("diff", "--stat"),
    "",
    "## 8. Dependency files"
]

for p in files:
    if p.name.lower() in {
        "requirements.txt", "pyproject.toml",
        "environment.yml", "setup.py"
    }:
        try:
            lines.append(f"\n### {p.relative_to(ROOT)}")
            lines.append(p.read_text(encoding="utf-8")[:10000])
        except Exception:
            pass

report = OUT / "project_audit.md"
report.write_text("\n".join(lines), encoding="utf-8")

archive = ROOT / "phyloresist_project_audit.zip"
with zipfile.ZipFile(archive, "w", zipfile.ZIP_DEFLATED) as z:
    z.write(report, "project_audit.md")

print("\nPROJECT AUDIT COMPLETED")
print(f"Files scanned: {len(files)}")
print(f"Audit report: {report}")
print(f"Upload this file: {archive}")
