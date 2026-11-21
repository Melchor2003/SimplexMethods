# Simplex Linear Programming Solver

An educational web app that solves linear programming problems with the simplex method and shows every step: tableaus, pivot elements, entering and leaving variables, ratios, and pivotal elimination.

Built for **Quantitative Methods** by **Group 3**, using plain HTML, CSS, and JavaScript with no frameworks and no build step.

## Pages

| Page | What it does |
| --- | --- |
| **Home** | Overview of the tool and its features |
| **Solver** | Enter an LP model and get a step-by-step solution |
| **About Team** | Team members, teacher, and project information |
| **Reviewer** | Simplex method summary, module downloads, and a feedback form |

## Features

- Maximize or Minimize objectives
- Constraints with ≤, ≥, and = (Big M method is applied automatically)
- Standard form shown before solving, including the LP Model line (for example `120x + 100y --> z - 120x - 100y = 0`)
- Non-negativity is always active and listed for every variable, including slack, surplus, and artificial ones
- Step-by-step tableaus with the pivot element, pivot column, and pivot row highlighted
- Ratio column, entering variable, and leaving variable for each iteration
- **Pivotal Elimination Steps** table for every iteration (new pivot row, other rows, and Z row, shown element by element)
- **Solution Check** that substitutes the answer back into every constraint, the non-negativity rule, and Z
- Input checking that highlights the field with the problem
- Fractions or decimals display toggle
- Detection of unbounded, infeasible, and alternative-optima cases
- Built-in examples grouped as Standard, Big M, and Special cases (includes the 120x + 100y class example)
- Light and dark theme, responsive layout

## Project structure

```
project-folder/
├── index.html     # all four pages (hash-routed: #/home, #/solver, #/team, #/reviewer)
├── style.css      # design tokens, layout, components, dark mode
├── script.js      # solver, routing, form handling, rendering, team, modules, feedback
├── README.md
├── images/        # team and teacher photos
└── assets/        # module files (PDF and PPT)
```

Keep all of these in the same folder.

## Run locally

Open `index.html` in any modern browser. No install or server is required.

## Deploy to Firebase Hosting

1. Install the CLI: `npm install -g firebase-tools`
2. Run `firebase login`, then `firebase init hosting` and choose `public` as the folder.
3. Copy `index.html`, `style.css`, `script.js`, the `images` folder, and the `assets` folder into `public/`.
4. Run `firebase deploy`.

## How to use the solver

1. Choose **Maximize** or **Minimize**.
2. Type the objective, for example `50x + 40y`. Variable names such as `x1`, `x2`, or `x`, `y`, `z` all work, and a coefficient of 1 can be left out.
3. Add constraints: left side (`2x + y`), operator, and right-hand side number (`100`).
4. Click **Solve LP Model**, or pick a model from **Load an example** first.

Example: Maximize `Z = 50x + 40y` subject to `2x + y ≤ 100`, `x + 2y ≤ 80` gives **Z = 2800** at x = 40, y = 20.

## Adding team and teacher photos

Photos load from the `images/` folder. The file names are the `photo` values in the `team` array (and `PROF`) in `script.js`. Do not type the extension in the code; jpg, jpeg, png, webp, and gif all work. If a photo is missing, the card shows the person's initials.

| Person | Photo file name |
| --- | --- |
| Rupha Alenna O. Lopega | `rupha` |
| Shiena M. Lucban | `Shiena` |
| Jethro G. Lucindo | `jethro` |
| Keith Ablow B. Mabanag | `keith` |
| Miguel Angelo R. Mercado | `miguel` |
| Andrea Paula C. Mopas | `paula` |
| Melchor M. Morada Jr. | `melchor` |
| Aljhen L. Nabus | `aljhen` |
| PTL. Norcelito Galvan (teacher) | `prof` |

File names are case-sensitive on most web hosts (including Firebase), so `Shiena.jpg` and `shiena.jpg` are different files. Match the capitalization exactly.

To change a name, role, location, quote, or photo name, edit the `team` array in `script.js`.

## Adding the module files (PDF / PPT)

Module downloads come from the `assets/` folder. The files and their descriptions are listed in the `MODULES` array in `script.js`.

| Module | File name |
| --- | --- |
| Module 1 – Linear Programming | `Module-1-Quantitative-Methods.pdf` |
| Module 2 – Simplex Method | `Module-2-QUANTITATIVE-METHODS2-1.pdf` |
| Module 3 – Two-Phase Method & M-Method | `Module-3-M-METHOD_AND_TP_Final.pdf` |
| Module 4 – Dual Method | `Module-4-The-Dual-Method.pdf` |
| Module 5 – Additional Materials | `qm-LP_090138.pptx` |

To add, rename, or remove a module, edit one line in `MODULES`. The file name there must match the file in `assets/` exactly.

## How the solver works

1. Constraints with a negative right-hand side are multiplied by −1 so every RHS is non-negative.
2. Slack variables are added for ≤, surplus and artificial variables for ≥, and artificial variables for =.
3. Artificial variables are penalized with Big M. The Z row stores each entry as `constant + coefficient·M`, so M is handled symbolically rather than as a large number.
4. Minimization is solved as Max Z′ = −Z, then converted back.
5. The entering column is the most negative Z-row value, and the leaving row is the smallest positive ratio.
6. The solver stops when no negative Z-row values remain. A positive artificial variable left in the basis means the problem is infeasible. A pivot column with no positive entries means it is unbounded.

## Notes

- Feedback on the Reviewer page is saved only in the visitor's own browser (`localStorage`). It is not sent to a server.
- Google Fonts (Inter and JetBrains Mono) are loaded from the web, with system fonts as a fallback.

## Team

| Member | Role | Location |
| --- | --- | --- |
| Rupha Alenna O. Lopega | Documentation | Gamay, Northern Samar |
| Shiena M. Lucban | UI/UX Designer | UEP Main Campus |
| Jethro G. Lucindo | Quality Assurance | Gamay, Northern Samar |
| Keith Ablow B. Mabanag | Back-End Developer | UEP Main Campus |
| Miguel Angelo R. Mercado | Mathematician | Catarman |
| Andrea Paula C. Mopas | Mathematician | UEP Main Campus |
| Melchor M. Morada Jr. | UI/UX Designer | Las Navas |
| Aljhen L. Nabus | Documentation | Pambujan |

**Teacher:** PTL. Norcelito Galvan, Quantitative Methods

© 2026 Quantitative Methods Group 3. All rights reserved.
