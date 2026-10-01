# Primary references and citation review

The following primary-source pages were checked on 2026-09-07. They support
the background methods, not this development version's tests or CoRE-MOF-COD
counts. Local implementation details are supported by [evidence.md](evidence.md).
Importable entries are provided in [references.bib](references.bib); abbreviated
author lists should be completed from the DOI records before submission.

| Key | Reference | Use in this draft |
|---|---|---|
| `coremof2025` | Zhao et al., “CoRE MOF DB: A curated experimental metal-organic framework database with machine-learned properties for integrated material-process screening,” *Matter* 8, 102140 (2025). [Publisher](https://www.sciencedirect.com/science/article/pii/S2590238525001833) | Database context; not a citation for the new 42,574-member release |
| `zeopp2012` | Willems et al., “Algorithms and tools for high-throughput geometry-based analysis of crystalline porous materials,” *Microporous and Mesoporous Materials* 149, 134–141 (2012). [Author-maintained method documentation](https://www.maciejharanczyk.info/Zeopp/docs.html) | Pore analysis; exact probe/settings remain release-specific |
| `rac2017` | Janet and Kulik, “Resolving Transition Metal Chemical Space: Feature Selection for Machine Learning and Structure–Property Relationships” (2017). [Publisher](https://doi.org/10.1021/acs.jpca.7b08750) | RAC method family; not a source for the project's exact 264-value schema |
| `mofid2019` | Bucior et al., “Identification Schemes for Metal–Organic Frameworks To Enable Rapid Search and Cheminformatics Analysis,” *Crystal Growth & Design* 19, 6682–6697 (2019). [Publisher](https://doi.org/10.1021/acs.cgd.9b01050) | Original identifier method; version-two extensions additionally require the applicable update/version citation |
| `crystalnets2022` | Zoubritzky and Coudert, “CrystalNets.jl: Identification of Crystal Topologies,” *SciPost Chemistry* 1, 005 (2022). [Publisher](https://doi.org/10.21468/SciPostChem.1.2.005) | Topology method; release fingerprint is separately defined |
| `sklearn2011` | Pedregosa et al., “Scikit-learn: Machine Learning in Python,” *JMLR* 12, 2825–2830 (2011). [Journal](https://www.jmlr.org/papers/v12/pedregosa11a.html) | Numerical implementation; exact package versions come from source/receipts |
| `umap2018` | McInnes, Healy, and Melville, “UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction,” arXiv:1802.03426 (2018; revised 2020). [Author manuscript](https://arxiv.org/abs/1802.03426) | Companion visualization; not a partitioning or duplicate criterion |

## Required author review before submission

- Complete citations for every actually executed checker, including the exact
  original MOFChecker, Chen–Manz, MOSAEC, MOFClassifier, and SETC-GAT methods.
  The [existing reference page](../docs/source/references.rst) is a starting
  list, not a substitute for checking version/method correspondence.
- If presenting newly calculated descriptors, document the actual molSimplify,
  Zeo++, CrystalNets, MOFid, and optional site-matching configurations from their
  receipts, with underlying software citations.
- Cite PACMAN, pretrained stability, and heat-capacity models only as used,
  with weight/version provenance. Ensemble dispersion is not automatically a
  calibrated prediction interval.
- The target-method supplement must document the actual RASPA version, force
  field, charge method, cutoffs, sampling, and endpoint definitions. A raw
  Widom weight ratio must not be renamed a finite-pressure mixture selectivity.
- Add the final software archive identifier and the appropriate approved
  dataset identifier. Neither has been invented for this documentation draft.
