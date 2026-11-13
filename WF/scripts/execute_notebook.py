"""
Executes a notebook in place with the kernel of the running Python environment (python3), from the notebook's folder,
and keeps the notebook's kernelspec metadata. Exit status 1 if a cell fails (the notebook is saved up to that cell).

Usage (e.g. in an sbatch job, sanity_utils environment):
    python WF/scripts/execute_notebook.py <notebook.ipynb>
"""

import os
import sys
import time

import nbformat
from nbclient import NotebookClient
from nbclient.exceptions import CellExecutionError

path = os.path.abspath(sys.argv[1])
nb = nbformat.read(path, as_version=4)
kernelspec = nb.metadata.get("kernelspec")
client = NotebookClient(nb, timeout=None, kernel_name="python3", resources={"metadata": {"path": os.path.dirname(path)}})
t0 = time.time()
try:
    client.execute()
    status = 0
except CellExecutionError as e:
    print(str(e)[:4000])
    status = 1
nb.metadata["kernelspec"] = kernelspec
nbformat.write(nb, path)
print(f"{os.path.basename(path)}: done in {time.time() - t0:.0f} s, status {status}")
sys.exit(status)
