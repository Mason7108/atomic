# Atomica — Interactive Atom Explorer

A self-contained, editable website featuring a searchable collection of all 118 chemical elements, simplified animated electron-shell diagrams, clickable protons/neutrons/electrons, a zoomable nucleus close-up, isotope selection, and learning notes.

## Try it

Double-click **`index.html`**. It works in a modern browser, without a server, internet connection, installs, or API keys.

## Edit it in VS Code

1. Open `index.html` in VS Code.
2. The design is in the `<style>` section. The interactivity is in the bottom `<script>` section.
3. Customize the colors, labels, animation speed, or add more explanations.
4. Save and refresh the browser. If you use VS Code Live Server, it will refresh automatically.

## Publish to GitHub Pages

1. Create a GitHub repository.
2. Add `index.html` to the repository root (the other files are optional).
3. In repository Settings → Pages, publish from the appropriate branch and the root folder.
4. Open the Pages URL GitHub provides.

## Features

- Type a **name**, **symbol**, or **atomic number** in search (e.g., `oxygen`, `O`, or `8`).
- Watch electrons move along an **illustrative Bohr-style orbital diagram** (not a physically accurate representation of quantum electron probability clouds).
- Inspect **pink protons**, **gold neutrons**, and **cyan electrons** and click them for descriptions.
- Switch to **Nucleus view** to see the individual nucleons on a larger scale.
- Drag on the canvas to adjust the illustrative rotation, scroll over it to zoom, or adjust the controls.
- Try selected known **isotopes**, such as Carbon-12 / Carbon-13 / Carbon-14 and Uranium-235 / Uranium-238.
- Pauses and reduces motion automatically if the device requests reduced motion.
- No dependencies: JS, CSS and element data are embedded in one file.

## Scientific notes

- Neutral atoms contain **Z protons and Z electrons**. This app depicts neutral atoms.
- For a selected isotope of mass number **A**, neutrons = **A − Z**. Displayed atomic weights are average/reference atomic weights and should not be treated as exact isotope mass numbers.
- The moving electrons and shell rings are a *learning diagram*, not actual electron orbits, sizes, distances, or speeds.
- Some electron configurations and properties of the heaviest elements are theoretical or simplified.

## Data

Element names, symbols, atomic weights and reference isotope numbers were generated from the **RDKit** periodic table. Electron shell counts were assembled from simplified orbital-filling rules with selected ground-state exceptions. The standalone app does not need RDKit to run. `elements.json` is supplied for reference, and `build_data.py` reproduces it if RDKit is installed.

