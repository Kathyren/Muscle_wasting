<h1>Muscle_wasting</h1>

<p>A bioinformatics pipeline for transcriptomic analysis of muscle wasting conditions. Integrates differential expression, transcription factor analysis, pathway enrichment, and network generation in a unified workflow.</p>

<h2>🚀 Features</h2>
<ul>
  <li><strong>run_DEA</strong> – Perform Differential Expression Analysis (DEA) and extract Differentially Dysregulated Sites (DDSs).</li>
  <li><strong>run_DDs</strong> – Identify Differentially Expressed Genes (DEGs).</li>
  <li><strong>run_TF</strong> – Predict and analyze transcription factors driving the DEGs.</li>
  <li><strong>run_pathway</strong> – Conduct pathway enrichment on DEGs or TF targets.</li>
  <li><strong>run_network_creation</strong> – Build regulatory networks based on TF–gene relationships.</li>
  <li><strong>main_pipeline</strong> – Execute the full workflow in sequence.</li>
</ul>

<h2>📁 Repository Structure</h2>
<pre><code>Muscle_wasting/
├── common_tools/          # Shared scripts and utilities
├── database_analysis/     # Database-driven analyses
├── data/                  # Raw/processed input data
├── network/               # Network construction modules
├── paper_mining/          # PubMed paper mining scripts
├── cytoscape/             # Cytoscape-specific tools
├── Databases/             # External data resources
├── ncbi/                  # NCBI query modules
├── Tests/                 # Unit and integration tests
├── Constants.py           # Global constants and config
├── main.py                # CLI entrypoint; defines main_pipeline
├── requirements.txt       # Python dependencies
└── DE_genes.ods           # Example DEG output file
</code></pre>

<h2>🛠️ Installation</h2>
<p>Install required packages:</p>
<pre><code>pip install -r requirements.txt</code></pre>

<p>Main dependencies:</p>
<ul>
  <li>pandas</li>
  <li>numpy</li>
  <li>scipy</li>
  <li>statsmodels (or DESeq2/edgeR via rpy2)</li>
  <li>networkx</li>
  <li>gseapy or enrichr</li>
  <li>py4cytoscape</li>
  <li>pytest</li>
</ul>

<h2>⚙️ Usage</h2>

<h3>As a library:</h3>
<pre><code>from main import run_DEA, run_DDs, run_TF, run_pathway, run_network_creation

dds = run_DEA(input_counts, sample_metadata)
degs = run_DDs(dds)
tfs = run_TF(degs)
pathways = run_pathway(degs)
network = run_network_creation(tfs, degs)
</code></pre>

<h3>Via the command-line:</h3>
<pre><code>python main.py \
  --input counts.csv \
  --metadata meta.csv \
  --output_dir results/
</code></pre>

<p>This runs the full <code>main_pipeline</code>, sequentially performing:</p>
<ol>
  <li>DEA & DDS detection</li>
  <li>DEG calling</li>
  <li>TF analysis</li>
  <li>Pathway enrichment</li>
  <li>Network construction</li>
</ol>

<p>Results are stored in <code>results/</code> as CSV and graph files.</p>

<h2>🧪 Testing</h2>
<p>Run tests from the repo root:</p>
<pre><code>pytest</code></pre>

<h2>📄 Output</h2>
<ul>
  <li>DDS list</li>
  <li>DEG list</li>
  <li>Predicted TFs and regulators</li>
  <li>Enriched pathway reports</li>
  <li>Network export formats (e.g. SIF, GraphML)</li>
</ul>

<h2>📚 Citation</h2>
<p>If you publish results using this pipeline, please cite the original <a href="https://github.com/GuerreroVazquez/Muscle_wasting">GuerreroVazquez/Muscle_wasting</a> repository.</p>

<h2>❓ Questions & Contributing</h2>
<p>For bugs, feature requests, or contributions, please open an issue or pull request.</p>

<hr />

<h3>Author</h3>
<p>
  Valeria Guerrero‑Vazquez<br />
  📧 <a href="mailto:GuerreroVazquez@gmail.com">GuerreroVazquez@gmail.com</a><br />
  GitHub: <a href="https://github.com/GuerreroVazquez">GuerreroVazquez</a>
</p>
