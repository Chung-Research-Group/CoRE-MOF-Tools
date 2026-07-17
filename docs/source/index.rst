.. title:: CoRE MOF Tools

.. raw:: html

   <section class="coremof-hero">
     <div class="coremof-hero__content">
       <span class="coremof-eyebrow">OPEN TOOLS FOR REPRODUCIBLE MOF RESEARCH</span>
       <h1>From MOF structures to research-ready insight.</h1>
       <p>
         A unified Python toolkit for accessing the CoRE MOF database, curating
         crystal structures, computing descriptors, and running validated
         property-prediction workflows.
       </p>
       <div class="coremof-actions">
         <a class="coremof-button coremof-button--primary" href="installation.html">Install CoRE MOF Tools</a>
         <a class="coremof-button coremof-button--secondary" href="quickstart.html">Explore the quick start</a>
       </div>
       <div class="coremof-meta" aria-label="Project highlights">
         <span>Python 3.9–3.11</span>
         <span>Command-line diagnostics</span>
         <span>Research-focused workflows</span>
       </div>
     </div>
     <div class="coremof-hero__brand" aria-hidden="true">
       <img src="_static/coremof-logo.png" alt="">
       <div class="coremof-lattice">
         <span></span><span></span><span></span><span></span><span></span>
       </div>
     </div>
   </section>

.. raw:: html

   <h2 id="workflows">One toolkit, four connected workflows</h2>

Move from source data to analysis without stitching together undocumented
scripts. Each workflow has explicit requirements, predictable outputs, and
troubleshooting guidance.

.. grid:: 1 2 2 4
   :gutter: 3
   :class-container: coremof-capability-grid

   .. grid-item-card:: Database access
      :link: quickstart
      :link-type: doc
      :class-card: coremof-card

      Query CoRE MOF records, inspect metadata, and retrieve associated
      adsorption data with clear dataset validation.

      +++
      **Start with:** ``information()``

   .. grid-item-card:: Curation & validation
      :link: features
      :link-type: doc
      :class-card: coremof-card

      Standardize CIFs, preserve ion information, and run structural checks
      before downstream calculations.

      +++
      **Outputs:** curated CIF + JSON report

   .. grid-item-card:: Geometry & descriptors
      :link: features
      :link-type: doc
      :class-card: coremof-card

      Calculate pore geometry, topology, open metal sites, and RAC descriptors
      through documented interfaces.

      +++
      **Integrates:** Zeo++, CrystalNets, molSimplify

   .. grid-item-card:: Property prediction
      :link: features
      :link-type: doc
      :class-card: coremof-card

      Apply packaged or repository-hosted models for stability, charge, and
      heat-capacity workflows.

      +++
      **Includes:** ensemble uncertainty summaries

.. raw:: html

   <h2 id="getting-started">Start with a healthy environment</h2>

Install in an isolated environment, then use the built-in diagnostic to see
which optional scientific components are available.

.. grid:: 1 1 2 2
   :gutter: 4
   :class-container: coremof-start-grid

   .. grid-item::

      .. code-block:: bash
         :caption: Terminal

         conda create -n coremof python=3.11
         conda activate coremof
         pip install CoREMOF-tools
         coremof doctor

   .. grid-item::

      .. raw:: html

         <div class="coremof-checklist">
           <div><strong>1</strong><span><b>Install the base package</b><small>Database lookup and diagnostics work independently of licensed tools.</small></span></div>
           <div><strong>2</strong><span><b>Run <code>coremof doctor</code></b><small>See feature availability and actionable setup guidance.</small></span></div>
           <div><strong>3</strong><span><b>Choose a workflow</b><small>Follow examples with documented units, outputs, and dependencies.</small></span></div>
         </div>

.. button-ref:: installation
   :ref-type: doc
   :color: primary
   :expand:

   Read the complete installation guide

.. raw:: html

   <h2 id="research-workflows">Built for real scientific workflows</h2>

.. grid:: 1 1 3 3
   :gutter: 3

   .. grid-item::
      :class: coremof-principle

      **Safer execution**

      External programs run without shell interpolation and use isolated
      temporary files for concurrent workloads.

   .. grid-item::
      :class: coremof-principle

      **Optional by design**

      Missing licensed or specialist software produces feature-specific,
      actionable errors instead of breaking unrelated imports.

   .. grid-item::
      :class: coremof-principle

      **Traceable results**

      Units, scientific parameters, model requirements, and generated files are
      documented where researchers need them.

.. note::

   CSD downloads and MOSAEC require a licensed CSD installation. Zeo++ and
   MOFid are external programs. Heat-capacity prediction requires ensemble
   files from the full GitHub repository. :doc:`installation` explains each
   optional component.

.. toctree::
   :maxdepth: 2
   :caption: User guide
   :hidden:

   installation
   quickstart
   features
   troubleshooting
   references

.. toctree::
   :maxdepth: 2
   :caption: API reference
   :hidden:

   CoREMOF

.. raw:: html

   <h2 id="project-resources">Project resources</h2>

.. grid:: 1 2 3 3
   :gutter: 2

   .. grid-item-card:: Source code
      :link: https://github.com/Chung-Research-Group/CoRE-MOF-Tools
      :link-type: url
      :class-card: coremof-resource

      Browse releases, implementation details, and contribution guidance.

   .. grid-item-card:: CoRE MOF database
      :link: https://mof-db.pusan.ac.kr/
      :link-type: url
      :class-card: coremof-resource

      Explore the associated database and web application.

   .. grid-item-card:: Report a problem
      :link: https://github.com/Chung-Research-Group/CoRE-MOF-Tools/issues
      :link-type: url
      :class-card: coremof-resource

      Submit a reproducible issue using the project template.

Indices: :ref:`genindex` · :ref:`modindex` · :ref:`search`
