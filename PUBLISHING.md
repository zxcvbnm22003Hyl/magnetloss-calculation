# GitHub publishing checklist

Before the first public release:

1. Choose the final repository name and GitHub organization/user.
2. Replace `APEX contributors` in `LICENSE`, `pyproject.toml`, and `CITATION.cff` with the intended copyright/author information if desired.
3. Check whether the project name, geometry, material table, and reference numerical values are approved for public disclosure.
4. Run:

   ```bash
   python -m pip install -e .[dev]
   pytest -q
   python scripts/validate_reference.py
   ```

5. Create the repository and push:

   ```bash
   git init
   git add .
   git commit -m "Initial open-source release v0.1.0"
   git branch -M main
   git remote add origin <repository-url>
   git push -u origin main
   ```

6. Create a GitHub release/tag `v0.1.0` and attach the wheel if you want a directly installable artifact.

The code is released under BSD-3-Clause in the current package. Change the license before publishing if your institution requires another license.
