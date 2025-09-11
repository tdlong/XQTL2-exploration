Production Directory Policy (scripts/production)

DO NOT add, modify, or delete files in this directory without explicit approval.

Enforcement in this repo:
1) Git pre-commit hook blocks changes to scripts/production/* by default.
2) CODEOWNERS requires @tdlong approval for any changes to scripts/production/*.
3) Local filesystem permissions are set read-only to prevent accidental edits.

To intentionally modify production (with approval):
- Temporarily bypass the pre-commit hook by exporting ALLOW_PROD_EDITS=1
  Example: ALLOW_PROD_EDITS=1 git commit -m "[prod-ok] approved change"
- After making changes, restore permissions using: chmod -R a-w scripts/production

Contact: @tdlong
