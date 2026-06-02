---
title: 'Panel Data Multiverse Lab: A Browser-Based Tool for Teaching Panel-Data Estimation and the Fragility of Empirical Results'
tags:
  - JavaScript
  - R
  - Python
  - econometrics
  - panel data
  - fixed effects
  - specification curve
  - multiverse analysis
  - reproducibility
  - teaching tool
authors:
  - name: Kevin Schoenholzer
    orcid: 0000-0001-9892-5869
    affiliation: 1
affiliations:
  - name: Institute of Communication and Public Policy, Università della Svizzera italiana (USI), Lugano, Switzerland
    index: 1
date: 26 May 2026
bibliography: paper.bib
---

# Summary

The **Panel Data Multiverse Lab** is a free, open-source, browser-based tool for teaching panel-data econometrics and the idea that an empirical estimate is one point in a space of defensible analytic choices. Working with a real longitudinal wage panel — 545 young men from the U.S. National Longitudinal Survey of Youth observed annually from 1980 to 1987 [@vellaverbeek1998], distributed in the `wooldridge` and `plm` R packages [@wooldridge2010pkg; @croissant2008] — students pick an estimand (the union wage premium, the marriage premium, or the return to schooling) and then build the entire **specification curve** [@simonsohn2020] implied by reasonable choices of estimator (pooled OLS, between, random effects, fixed effects, first differences, correlated random effects, and two-way fixed effects), control set, sample, and standard-error family. The application performs all computation client-side in dependency-free JavaScript, requires no installation, no server, and no programming, and runs on any device with a web browser. Every estimator is implemented from textbook formulas and **independently verified against R's `plm` and Python's `linearmodels`** across all 1,440 specifications of the default multiverse. The tool is accompanied by a live "under-the-hood" panel that displays the equation each estimator is fitting and updates it as choices change, a one-click translation of any specification into runnable R, Stata, and Python code, exportable figures and reports, a guided lesson, and self-check questions.

# Statement of need

Empirical economics has reorganized itself around research design [@angrist2010credibility], yet introductory and intermediate econometrics is still largely taught as a sequence of single point estimates. Students too often infer that the data deliver one correct number, when practising researchers know an estimate is the endpoint of a long chain of analytic decisions — which estimator, which controls, which sample, which standard errors — many of them defensible yet consequential. Naming this problem has a long history, from Leamer's [-@leamer1983] call to "take the con out of econometrics" to modern evidence of specification searching [@brodeur2016] and the tools developed to confront it: the garden of forking paths [@gelman2013], multiverse analysis [@steegen2016], and specification-curve analysis [@simonsohn2020].

Teaching this idea is hard because constructing a specification curve normally requires programming fluency that beginning students do not yet have, so the lesson is deferred until after the habit of the single estimate has formed. Existing resources do not close the gap: specification-curve packages such as `specr` and panel tools such as `plm` [@croissant2008] or Stata's `xtreg` are code-first and presuppose the skills the lesson is meant to build, while most interactive statistics applets rely on simulated data that students can dismiss as artificial. Panel data are an especially sharp setting for the lesson, because the choice of estimator encodes a well-understood identification trade-off [@mundlak1978; @wooldridge2010]: fixed effects remove all time-invariant confounding but cannot identify a regressor that never changes within a unit. The Panel Data Multiverse Lab fills this gap with a no-code, real-data, verified environment in which students build the whole universe of estimates before committing to any one, and read identification and fragility directly off the result.

# Functionality

The application is organized into linked views:

- **Data overview** presents the panel as people (individual wage trajectories, distributions, a within- versus between-person variance decomposition, and a table of "switchers" that makes concrete which observations identify a fixed-effects estimate).
- **Specification curve** enumerates every combination of the selected analytic choices (typically a few hundred), sorts the estimates, and plots them above a dashboard of the choices behind each, with summary statistics (median, range, share significant, share sign-flipped, number unidentifiable) and an optional randomization test of the whole curve.
- **Single-spec deep-dive** opens any one specification with a full regression table, all estimators side by side, a Hausman test, residual diagnostics, a cluster bootstrap, and a live equation panel that shows the model each estimator fits and the standard-error formula in use.
- **Export & code** translates the current specification into runnable R (`plm`), Stata, and Python (`linearmodels`) code and exports the multiverse as CSV, the figures as SVG/PNG, and a self-contained HTML report.
- **Methodology** documents the data, estimators, inference, verification, a glossary, self-check questions, and citation information.

A scripted guided lesson, per-tab learning objectives, a predict-then-reveal interaction, and shareable permalinks that encode a full specification in the URL support both self-study and instructor-led use.

# Verification

Because a teaching tool that produced subtly wrong numbers would mislead the students it aims to help, the numerical core is held to a research standard of reproducibility. A self-test reproduces canonical `plm` benchmarks, and every estimate in the default multiverse — all 1,440 specifications across three estimands — is re-fit independently in two languages, R (`plm`) and Python (`linearmodels`). The pooled-OLS, between, fixed-effects and first-difference coefficients agree to better than $1\times10^{-6}$ (R) and $5\times10^{-7}$ (Python), and random effects to $3\times10^{-4}$; the Mundlak (correlated random effects) coefficient reproduces the fixed-effects estimate exactly, as theory requires. The data, source, and all three verification scripts are openly available.

# Target audience

The tool is designed for intermediate-undergraduate through graduate courses in econometrics, applied microeconomics, and quantitative social science, at the point where panel estimators are introduced. It assumes familiarity with multiple regression but no programming, which lets the conceptual lesson precede the acquisition of software skills, and it can be used as a lecture demonstration, a guided lab, or a bridge to reproducible coding.

# Acknowledgements

The wage-panel extract derives from @vellaverbeek1998 and is redistributed through the `wooldridge` [@wooldridge2010pkg] and `plm` [@croissant2008] R packages. The tool was developed by the author with AI-assisted coding.

# References
