# Deploying this site with GitHub Pages

This repository includes a simple GitHub Pages site scaffold in the `docs/` folder. Follow these steps to get it live:

1. Add your real files:
   - Put `writeup.pdf`, `slides.pdf`, and `software.tar.gz` into `docs/assets/`.
   - Replace the placeholder portraits in `docs/index.html` with images in `docs/images/` (create `docs/images/` and commit images there), then update the `src` attributes.

2. Commit and push the changes to GitHub:

```bash
git add docs/ && git commit -m "Add website scaffold" && git push origin main
```

3. Enable GitHub Pages (docs folder) in the repository settings:
   - On GitHub, go to Settings -> Pages
   - Under "Source", choose: Branch: `main`, Folder: `/docs`
   - Save. The page will be published at `https://<your-username>.github.io/<repo-name>/`.

Alternative: publish from a `gh-pages` branch. If you prefer that workflow, create a branch and push built site there.

Notes:
- If you want HTTPS or custom domain, configure those in Pages settings.
- To preview locally, open `docs/index.html` in a browser. For more advanced previews use a local static server, e.g.: `python -m http.server --directory docs 8000` and visit `http://localhost:8000`.
