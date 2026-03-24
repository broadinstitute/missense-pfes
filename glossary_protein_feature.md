# Glossary - protein features

## **Physicochemical properties of amino acids (`Physicochemical`)**

Amino acids are grouped based on the chemical nature of their side chains, which influences protein structure, stability, and interactions.

### **Reference amino acid classes**

* **Aliphatic** – Nonpolar and hydrophobic amino acids with straight or branched hydrocarbon side chains.
  * Alanine (Ala, A), Isoleucine (Ile, I), Leucine (Leu, L), Methionine (Met, M), Valine (Val, V).
* **Aromatic** – Contain an aromatic ring; often involved in stacking interactions.
  * Phenylalanine (Phe, F), Tryptophan (Trp, W), Tyrosine (Tyr, Y)
* **Polar/Neutral** – Polar but uncharged side chains; can form hydrogen bonds.
  * Asparagine (Asn, N), Glutamine (Gln, Q), Serine (Ser, S), Threonine (Thr, T)
* **Positively-Charged** – Carry a positive charge at physiological pH; often interact with DNA or negatively charged groups
  * Arginine (Arg, R), Histidine (His, H), Lysine (Lys, K)
* **Negatively-Charged** – Carry a negative charge at physiological pH; participate in salt bridges and electrostatic interactions.
  * Aspartic acid (Asp, D), Glutamic acid (Glu, E)
* **Special** – Have unique structural or chemical roles that set them apart.
  * Proline (Pro, P) – Causes rigid kinks due to its cyclic structure.
  * Glycine (Gly, G) – Smallest residue, provides flexibility.
  * Cysteine (Cys, C) – Contains a thiol group, can form disulfide bonds.

### **Change in amino acid class**

Indicates the shift in physicochemical property classes caused by amino acid substitution (e.g., missense mutations):

* S441N - **Polar/Neutral>Polar/Neutral**
* G124D - **Special>Negatively-Charged**

### **Grantham's distance -- Severity of physicochemical perperty shift**

This classification is based on **Grantham's distance (*D*)**, a numerical measure of the biochemical difference between two amino acids. It considers **composition, polarity, and molecular volume** to quantify how drastic a substitution is. A higher *D* indicates a more disruptive and significant change in side chain chemistry. 

* **Mild** (*D* < 50) – Minimal biochemical change
* **Moderate** (50 ≤ *D* < 100) – Noticeable but not drastic change; slight shift in polarity or volume
* **Substantial** (100 ≤ *D* < 150) – Significant difference in chemical properties
* **Severe** (*D* ≥ 150) – Extreme shift with large biochemical mismatch

## **Functional Sites (`Function`)**

Annotations from UniProtKB describing key protein functional elements. See https://www.uniprot.org/help/function_section for more information. 

* **Active site** – Catalytic site of an enzyme.
* **Binding site** – Location where ligands (ions, cofactors, substrates) bind.
* **Site** – Single amino acid sites of functional interests; frequently used for annotating cleavage sites, Inhibitory sites for proteases, breakpoint sites for fusion proteins due to chromosomal rearrangement.
* **DNA binding** – Region interacting with DNA.
* **Zinc finger** – Structural motif coordinating one or more zinc ions.

## **Domains & Regions (`Domain`)**

Domain and region annotations are collected from UniProtKB to describe structural and functional segments of proteins.

* **Domain** – The position and type of each modular protein domain
* **Topological domain** – Location of non-membrane regions of membrane-spanning proteins (cytoplasmic or lumena sides)
* **Transmembrane** – Regions spanning the lipid bilayer membrane.
* **Intramembrane** – Regions embedded within the membrane but not spanning it entirely.
* **Compositional bias** – Region of compositional bias in the protein
* **Repeat** – Repeated sequence motifs or repeated domains
* **Motif** – Short (up to 20 amino acids) sequence motif of biological interest
* **Coiled coil** – Positions of regions of coiled coil within the protein; Structural motif facilitating oligomerization
* **Signal peptide** – Sequence targeting proteins to the secretory pathway or periplasmic space
* **Transit peptide** – Extent of a transit peptide for organelle targeting
* **Propeptide** – Part of a protein that is cleaved during maturation or activation
* **Peptide** – Extent of an active peptide in the mature protein
* **Chain** – Extent of a polypeptide chain in the mature protein

UniProtKB’s “Region” annotation covers areas of interest not captured by other categories. Because these can vary in nature, we subdivide them as follows:

* **Region/Disordered** – Regions annotated as ‘disordered,’ indicating intrinsically disordered regions (IDRs).
* **Region/Interaction** – Regions annotated for interaction with other molecules or proteins.
* **Region/Others** – Miscellaneous regions with other specific annotations.

See https://www.uniprot.org/help/family_and_domains_section for more information. 

## **Post-translational modifications (`Modification`)**

Annotations from UniProtKB and PhosphoSitePlus database describing post-translational modifications. This is a category of protein features that includes chemical modifications that occur after translation, which are essential for the proper folding, localization, activity and stability of the mature protein.

From UniProtKB (https://www.uniprot.org/help/ptm_processing_section):

* **Lipidation** – Covalently attached lipid group(s).
* **Glycosylation** – Covalently attached glycan group(s).
* **Crosslinks** – Residues participating in covalent linkage(s) between proteins.
* **Modified residue** – Modified residues excluding lipids, glycans and protein crosslinks.
* ~~**Disulfide bond (UniProtKB)** – Both experimentally determined and predicted disulfide bonds are annotated. it may refer to interchain or intrachain context depending on annotation (practically, it's very likely to be intrachain; > 99% intra, < 1% inter). <span style="color: red;">(*Aug 15th - merged to intra- and inter-molecular disulfide bond*)</span>~~ 

From PhosphositePlus (https://www.phosphosite.org/):

* **Acetylation** – Addition of an acetyl group, often on lysine residues; regulates protein stability and gene expression.
* **Methylation** – Addition of a methyl group, commonly on lysine or arginine; influences transcription and chromatin structure.
* **Phosphorylation** – Addition of a phosphate group (usually on Ser, Thr, or Tyr); key regulator of signaling pathways.
* **SUMOylation** – Attachment of SUMO (Small Ubiquitin-like Modifier) proteins; affects localization, stability, and activity.
* **Ubiquitination** – Covalent addition of ubiquitin; often tags proteins for degradation via the proteasome.
* **O-GalNAc/GlcNAc** – Specific O-linked glycosylation involving N-acetylgalactosamine or N-acetylglucosamine; modulates protein function and signaling.

## **3D Structural Features (`Structure`)**

Structural features are collected and calculated based on AlphaFold2 (secondary structures, solvent accessilbility, intramolecular interactions) and available PDB structures (intramolecular interactions).

### **DSSP Secondary Structure**

The classification of local secondary structure elements in a protein, as defined by the **DSSP (Define Secondary Structure of Proteins) algorithm**. DSSP assigns each amino acid residue to a structural category based on hydrogen bonding patterns and geometrical features derived from 3D protein structures. This is done on AlphaFold2 structures.

* **H** – Alpha-helix
* **G** – 3₁₀-helix
* **P** – Polyproline helix
* **I** – π-helix
* **B** – Isolated beta-bridge
* **E** – Extended beta-strand (in beta-sheet)
* **T** – Turn
* **S** – Bend
* **C** – Coil/undefined

### **Relative Solvent Accessibility (RSA)**

The **RSA** quantifies how exposed an amino acid residue is to solvent within a given structure. DSSP algorithm calculates the **absolute solvent accessible surface area** (ASA, measured in Å) of a residue *X* for a given protein structure. RSA is defined as ASA normalized against its maximum possible exposure, which is the ASA of the residue X in a reference tripeptide (Gly–*X*–Gly) configuration. These reference values were collected from the literature (*M. Z. Tien et al.*, *PLoS ONE* 8, e80635 (2013)).

Based on RSA values, residues were categorized into five exposure levels:

* **Core** – RSA < 5%
* **Buried** – 5% ≤ RSA < 25%
* **Medium-buried** – 25% ≤ RSA < 50%
* **Medium-exposed** – 50% ≤ RSA < 75%
* **Exposed** – RSA ≥ 75%

### **Structural Confidence (pLDDT)**

The **predicted Local Distance Difference Test (pLDDT)** is a per-residue confidence score from AlphaFold2, ranging from 0 to 100. It estimates the local accuracy of predicted structures without requiring alignment to experimental data. Regions with low pLDDT may reflect intrinsic disorder or insufficient evolutionary constraints.

* **Very high (pLDDT > 90)** – High confidence in backbone and side chains.
* **Confident (70 < pLDDT ≤ 90)** – Reliable backbone; side chains may be inaccurate.
* **Low (50 < pLDDT ≤ 70)** – Uncertain local structure.
* **Very low (pLDDT ≤ 50)** – Likely disordered or flexible regions.

### **Intramolecular Interactions**

Intra-protein interactions are detected from AlphaFold2 and monomeric PDB structures. Interactions present within a protein play critical roles in stabilizing its 3-dimensional structures.

* **Intramolecular hydrogen bond** – Hydrogen bond formed within a single protein chain.
* **Intramolecular salt bridge** – Electrostatic interaction between oppositely charged residues in the same protein chain.
* **Intramolecular disulfide bond** – Covalent bond between cysteine residues within a single protein chain. 
* **Intramolecular non-bonded interaction** – Hydrophobic or van der Waals contacts within a single protein structure. These are generally weak and can be transient. 

For interactions derived from AlphaFold2, we applied an additional filter using **Predicted Aligned Error (PAE)**, which is provided along with the AlphaFold2 prediction. PAE measures the confidence in the relative positions of two residues within the predicted structure and is defined as the expected positional error measured in Angstrom (Å). A high PAE score indicates the low confidence in their relative distance. As such, we excluded potential interactions between residues if their PAE is larger than 8Å, as there is a higher chance that such contacts are due to uncertainty in the model rather than true structural proximity. 

Please note that these interactions are inferred from distance-based criteria (see https://github.com/broadinstitute/g2p-bis/ for detail) applied to *static* protein structure and do not necessarily imply functional or energetic significance. The actual interaction pattern may vary due to protein flexibility, ligand binding, or dynamic changes in structure under physiological conditions.

Intramolecular disulfide bond also includes 'Disulfide bond' sequence annotation found in UniProtKB when annotated as intramolecular disulfide bonds. 

## **Protein-Protein Interactions (`PPI`)**

From experimentally resolved protein complex structure, interactions between two proteins found at their interface are detected. 

* **Intermolecular hydrogen bond** – Hydrogen bond formed between residues of different protein chains.
* **Intermolecular salt bridge** – Electrostatic interaction between oppositely charged residues on different protein chains.
* **Intermolecular disulfide bond** – Covalent bond formed between cysteines from separate protein chains.
* **Intermolecular non-bonded interaction** – Includes hydrophobic or van der Waals contacts between residues across protein interfaces.

Please note that these interactions are inferred from distance-based criteria (see https://github.com/broadinstitute/g2p-bis/ for detail) applied to static structural data and do not necessarily imply functional or energetic significance. The actual interaction pattern may vary due to protein flexibility, binding events, or dynamic changes in structure under physiological conditions.

Intermolecular disulfide bond also includes 'Disulfide bond' sequence annotation found in UniProtKB when tagged with 'Interchain'. 
