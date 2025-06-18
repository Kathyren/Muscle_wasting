<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
</head>
<body>
  <h1>miRKat: microRNA Toolkit for Muscle Wasting and Aging</h1>

  <p><strong>miRKat</strong> is a modular toolkit for exploring, evaluating, and prioritizing microRNAs (miRNAs) in muscle-related contexts such as wasting, regeneration, and aging. The pipeline integrates transcriptomic data, transcriptional networks, pathway analysis, and miRNA impact evaluation to help identify regulatory hotspots.</p>

  <h2>🧭 Overview</h2>
  miRKAt is build with three main parts:
Pre-processing (mirkitten)
miRKat Network
miRKat Scoring

<h2>Pre-processing (miRKitten)</h2>
  <ol>
    <li><strong>Pre-processing (mirkitten/)</strong><br>
        Generate networks by combining DE analysis, TF inference, and pathway enrichment.
        <ul>
          <li><code>run_DDS.py</code> – Combine DEA results</li>
          <li><code>run_DEG.py</code> – Identify differentially expressed genes</li>
          <li><code>run_TF.py</code> – Rank transcription factors</li>
          <li><code>run_pathway.py</code> – Perform pathway enrichment</li>
          <li><code>run_network_creation.py</code> – Generate gene-TF-microRNA network</li>
        </ul>
        These scripts can be run individually or via the full <code>network_creation</code> script.
    </li>

    <li><strong>Network Evaluation (network/main_pipeline.py)</strong><br>
        Select relevant subnetworks by defining criteria in a YAML config or via command-line options. Includes filters for pathways, TFs, and gene-level scores. 
    </li>

    <li><strong>miRNA Scoring (mirna_scoring/)</strong><br>
        Analyse the influence of miRNAs using a combination of graph structure, differential expression, pathway associations, and clustering. The scoring functions are callable via Python, with usage examples in:
        <br><a href="https://github.com/GuerreroVazquez/Muscle_wasting/blob/GuerreroVazquez-master/network/mirnas_influence_00.ipynb" target="_blank">mirnas_influence_00.ipynb</a>
    </li>
  </ol>

  
  <h2>⚙️ create_network Parameters</h2>
<p>The <code>create_network</code> function (used by <code>main_pipeline.py</code>) accepts a number of customizable parameters either via command-line arguments or a YAML configuration file. These define how the network is built, filtered, and enriched.</p>

<table border="1" cellpadding="6" cellspacing="0">
  <thead>
    <tr>
      <th>Parameter</th>
      <th>Description</th>
      <th>Example</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td><code>dds_threshold</code></td>
      <td>Minimum expression fold-change (log2) for gene filtering.</td>
      <td><code>2</code></td>
    </tr>
    <tr>
      <td><code>top_n_deg</code></td>
      <td>Number of top differentially expressed genes to retain.</td>
      <td><code>300</code></td>
    </tr>
    <tr>
      <td><code>tf_score_cutoff</code></td>
      <td>Minimum score for keeping transcription factors (used after TF ranking).</td>
      <td><code>0.3</code></td>
    </tr>
    <tr>
      <td><code>pathway_keywords</code></td>
      <td>List of keywords used to select enriched biological pathways.</td>
      <td><code>["MUSCLE", "MITO", "APOPTOSIS"]</code></td>
    </tr>
    <tr>
      <td><code>include_mirna</code></td>
      <td>Whether to include microRNAs in the network.</td>
      <td><code>true</code></td>
    </tr>
    <tr>
      <td><code>cutoff_score</code></td>
      <td>Filter for edge confidence or importance in final network.</td>
      <td><code>0.95</code></td>
    </tr>
    <tr>
      <td><code>input_deg_file</code></td>
      <td>Path to a CSV file containing pre-computed DEGs (instead of running DEA).</td>
      <td><code>"data/deg_results.csv"</code></td>
    </tr>
    <tr>
      <td><code>expression_profiles</code></td>
      <td>Path to expression matrix (tissue/cell-specific filtering).</td>
      <td><code>"data/protein_atlas_tissue.csv"</code></td>
    </tr>
    <tr>
      <td><code>cell_type_expression</code></td>
      <td>Cell-type expression data (optional).</td>
      <td><code>"data/muscle_cell_atlas.csv"</code></td>
    </tr>
    <tr>
      <td><code>save_path</code></td>
      <td>Output path to save the resulting network object (.pkl).</td>
      <td><code>"Networks_pkl/network_case5.pkl"</code></td>
    </tr>
    <tr>
      <td><code>usecase_name</code></td>
      <td>Label to group outputs from the same run.</td>
      <td><code>"UseCase_5"</code></td>
    </tr>
  </tbody>
</table>

<p>Parameters can be passed in a YAML config like:</p>

<pre><code>dds_threshold: 2
top_n_deg: 300
tf_score_cutoff: 0.3
pathway_keywords:
  - MUSCLE
  - MITO
  - APOPTOSIS
include_mirna: true
cutoff_score: 0.95
save_path: network/Networks_pkl/case5_network.pkl
usecase_name: UseCase_5
</code></pre>

<p>To get the dynamic DE threshold for all the comparisons using stat:</p>
<pre><code>python mirkitten/run_network_creation.py --dds_files UseCases/UseCase0/dds_files.yml --save_name UseCases/UseCase0/usecase0.cyjs --only_DE True</code></pre>

<p>To get the dynamic DE threshold for all the comparisons using log2FoldChange:</p>

<pre><code>python mirkitten/run_network_creation.py --dds_files UseCases/UseCase0/dds_files.yml --save_name UseCases/UseCase0/usecase0.cyjs --only_DE True --interest log2FoldChange</code></pre>

<h2>mirKat Network</h2>

<h3>Description</h3>
<p>
This module generates a biologically relevant network connecting <strong>miRNAs</strong>, <strong>mRNAs</strong>, and <strong>transcription factors (TFs)</strong> by integrating curated interactions, differential expression data, pathway associations, tissue specificity, and scoring strategies. It builds on a Cytoscape-compatible <code>.cyjs</code> base file, annotates nodes with relevant biological metadata, scores interactions based on relevance to the experimental context, and filters the network based on thresholds set by the user.
</p>

<h3>Parameters</h3>
<table border="1" cellpadding="6" cellspacing="0">
  <tr>
    <th>Parameter</th>
    <th>Type</th>
    <th>Description</th>
  </tr>
  <tr>
    <td><code>--input_cyjs</code></td>
    <td>str</td>
    <td>Path to the base network in <code>.cyjs</code> format.</td>
  </tr>
  <tr>
    <td><code>--tissue</code></td>
    <td>str</td>
    <td>Target tissue or cell type for annotation (e.g., "Vastus Lateralis").</td>
  </tr>
  <tr>
    <td><code>--pathway_db</code></td>
    <td>str</td>
    <td>Pathway database to use for scoring relevance (e.g., KEGG, Reactome).</td>
  </tr>
  <tr>
    <td><code>--expression_file</code></td>
    <td>str</td>
    <td>CSV file with differential expression scores (gene or miRNA-level).</td>
  </tr>
  <tr>
    <td><code>--score_weights</code></td>
    <td>dict</td>
    <td>User-defined weights for scoring features (e.g., pathway match, DEG status).</td>
  </tr>
  <tr>
    <td><code>--tf_database</code></td>
    <td>str</td>
    <td>Optional TF annotation file (if not in SQL backend).</td>
  </tr>
  <tr>
    <td><code>--filter_threshold</code></td>
    <td>float</td>
    <td>Cutoff value to retain nodes based on their final score.</td>
  </tr>
  <tr>
    <td><code>--output_name</code></td>
    <td>str</td>
    <td>Prefix for output files (e.g., "VL_miRNet").</td>
  </tr>
</table>

<h3>Usage Example</h3>
<pre><code>python main_pipeline.py \
    --input_cyjs base_network.cyjs \
    --tissue "Vastus Lateralis" \
    --pathway_db reactome.db \
    --expression_file VL_DEG_scores.csv \
    --score_weights '{"deg": 1.5, "pathway": 1.0, "tissue": 1.2}' \
    --tf_database tf_data.csv \
    --filter_threshold 2.0 \
    --output_name VL_miRNet
</code></pre>

<h3>Output Description</h3>
<ul>
  <li><code>VL_miRNet_filtered.cyjs</code>: Cytoscape-ready filtered and weighted network.</li>
  <li><code>VL_miRNet_scores.csv</code>: CSV file with node and edge scores and metadata annotations.</li>
  <li><code>VL_miRNet_report.txt</code>: Summary log with counts of nodes, edges, filtering steps, and key statistics.</li>
  <li>Optional visualizations (if enabled): subnetwork plots, score distributions, or TF–miRNA–gene hubs.</li>
</ul>

<h2> miRKAt Scoring</h2>
This is the final part, once you have all your personalized data an your sall network is time to see the microRNAs that impact the most your interests. You can do random walks to calculate the impact of the microRNAs on the network and find clusters to helop you to see redundancy on the mirnas. 
The results can be generated on html for you to see them separetly. An example of use can be seen on https://github.com/GuerreroVazquez/Muscle_wasting/blob/GuerreroVazquez-master/mirna_scoring/notebooks/Evaluate_mirs/UseCase_n.ipynb



  <h2>📊 Optional DEA</h2>
  <p>While DEA (Differential Expression Analysis) based on age groups (young, middle-aged, old) can be performed with <code>run_DEA.py</code>, it is <em>not</em> a required step for using miRKat. It is included to support aging-related use cases, such as muscle aging clocks.</p>

  <h2>🔍 Customization</h2>
  <p>miRKat supports personalisation through user-provided expression data for tissues and cell types:</p>
  <ul>
    <li>Tissue data: <a href="https://www.proteinatlas.org/humanproteome/tissue" target="_blank">Human Protein Atlas</a></li>
    <li>Muscle cell types: <a href="https://www.muscleageingcellatlas.org/" target="_blank">Human Muscle Cell Atlas</a></li>
  </ul>

  <h2>📁 Repository Structure</h2>
  <ul>
    <li><code>mirkitten/</code>: Pre-processing scripts</li>
    <li><code>network/</code>: Evaluation of subnetworks and filtering</li>
    <li><code>mirna_scoring/</code>: Analysis and scoring of miRNA influence</li>
    <li><code>UseCases/</code>: YAML configs and examples for running workflows</li>
  </ul>

  <h2>📌 Example Usage for miRNA Scoring</h2>
  <pre><code># Load and evaluate a network
import mirna_scoring.mirna_impact as mis
from network.network_processing import load_graph

network = load_graph("network/Networks_pkl/complete_n_tf_mirnas__UseCase5_cutoff_0.95.pkl")
evaluation = mis.mirna_evaluation(mis.mirna_network(network), session_name='UseCase_5')

# Score miRNAs
evaluation_df = evaluation.score(steps=10, sample_size=10, dds_threshold=2, 
                                 pathway_keywords=["MUSCLE", "MITO", "APOPTOSIS"])

# Cluster and plot
evaluation.cluster_mirnas(n_clusters=len(evaluation.get_all_mirnas())//3)
</code></pre>

  <p>See the notebook for full visualization steps using UMAP, pathway enrichment, and final reporting.</p>

  <h2>📜 License</h2>
  <p>MIT License</p>
</body>
</html>
