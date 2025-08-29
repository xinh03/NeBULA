# NeBULA
A novel drug design tool via bioisosteric replacement

![image](https://github.com/xinh03/NeBULA/blob/main/Figure_github.png)

<h2>Requirements</h2>
<pre><code>conda create -n nebula python=3.9
pip install numpy==1.21.6
pip install rdkit==2022.3.5
pip install pandas</code></pre>
<p>Note: you can use -- pip install "numpy<2.0" -- for more compatible version</p>

<h2>Folder format</h2>
<pre><code>--nebula
    --reaction.py
    --Reaction_Fsp3-rich.csv
    --README.md</code></pre>

<h2>How to use</h2>
<pre><code>
conda activate nebula
cd nebula
</code></pre>
<p>Note: you can use -- unzip Reaction_Fsp3-rich.zip -- to unzip necessary file</p>

<h3>1. Pass the parameter (in non-interactive mode)</h3>
<p>Filename <em>reaction.py</em>, which can be called in the terminal like this:</p>
<pre><code>python reaction.py --smiles "c1(c(COc2nc([C@H]3CC[N@@H+](Cc4nc5ccc(C(=O)[O-])cc5n4C[C@@H]4CCO4)CC3)ccc2)ccc(C#N)c1)F"</code></pre>
<p>"c1(c(COc2nc([C@H]3CC[N@@H+](Cc4nc5ccc(C(=O)[O-])cc5n4C[C@@H]4CCO4)CC3)ccc2)ccc(C#N)c1)F" means the SMILES you want to replace.</p>

<h3>2. Without passing parameters (interactive mode)</h3>
<p>If executed directly in the terminal,</p>
<pre><code>python reaction.py</code></pre>
<p>the program prompts for input:</p>
<pre><code>Enter SMILES:</code></pre>
<p>Input the appropriate content and press Enter to execute, the program will be based on the input to generate a "current time folder" contains: </p>
<pre><code>--%Y%m%d_%H%M%S
    --reaction_products.csv</code></pre>

