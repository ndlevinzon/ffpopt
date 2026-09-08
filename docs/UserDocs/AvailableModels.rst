Available Models
================

Many of the ffpopt scripts have a ``--model=str`` option.
The models that are available depend on how you installed
ffpopt because the backend libraries are not mutually
compatible with each other.

- sander
     - Description: MM calculation
     - Link: https://ambermd.org/
     - Note: This interface requires parm7/rst7 files

- uff
     - Description: MM calculation using the UFF force field
     - Note: Implemented within RDKit
     - Note: Prefer if metal ions are present

- mmff94
     - Description: MM calculation using the MMFF94 force field
     - Note: Implemented within RDKit
     - Note: Prefer if drug like organic molecule

- dftb2
     - Description: Second order self-consistent charge
       density functional tight binding. MIO-1-1 parameters.
     - Link: https://github.com/dftbparams/mio
     - Elements: H,C,N,O,P,S,Si,Ag,Ga
     - Note: This is run through pysander; therefore it
       requires parm7/rst7 files

- dftb3
     - Description: Full 3rd-order SCC-DFTB with O3B parameters.
     - Link: https://github.com/dftbparams/3ob
     - Elements: H,C,N,O,F,K,P,S,F,Na,Mg,Cl,Ca,I,Br,Zn
     - Note: This is run through pysander; therefore it
       requires parm7/rst7 files

- qdpi2
     - Description: xtb+delta MLP model based on DeepPot-SE
     - Link: https://www.doi.org/10.1021/acs.jpcb.4c01466
     - Elements: H,C,N,O,F,Na,P,S,Cl,K,Br,I
     
- xtb
     - Description: The GFN2-xTB semiempirical model
     - Link: https://github.com/grimme-lab/xtb
     - Link: https://pubs.acs.org/doi/10.1021/acs.jctc.8b01176
     - Elements: 1-84 (H-to-Po)
- mace-off23b_medium
     - Description: MACE-OFF23b_medium.model GNN
     - This is like mace-off23_medium; however, it increases the
       graph cutoff radius from 5A to 6A. This was found to be
       crucial for recovering certain condensed phase
       properties.
     - Link: https://github.com/ACEsuit/mace-off
     - Elements: H, C, N, O, F, P, S, Cl, Br, I

- mace-off23_small
     - Description: MACE-OFF23_small.model GNN
     - Link: https://github.com/ACEsuit/mace-off
     - Elements: H, C, N, O, F, P, S, Cl, Br, I

- mace-off23_large
     - Description: MACE-OFF23_large.model GNN
     - Link: https://github.com/ACEsuit/mace-off
     - Elements: H, C, N, O, F, P, S, Cl, Br, I

- mace-off24_medium
     - Description: MACE-OFF23_large.model GNN
     - Link: https://github.com/ACEsuit/mace-off
     - Elements: H, C, N, O, F, P, S, Cl, Br, I

- aimnet2
     - Description: Trained against wB97M-D3
     - Link: https://github.com/isayevlab/aimnetcentral
     - Elements: H, B, C, N, O, F, Si, P, S, Cl, As, Se, Br, I

- aimnet2_b973c
     - Description: Trained against B97-3c
     - Link: https://github.com/isayevlab/aimnetcentral
     - Elements: H, B, C, N, O, F, Si, P, S, Cl, As, Se, Br, I

- aimnet2_2025
     - Description: Trained against B97-3c + improved intermolecular interactions
     - Link: https://github.com/isayevlab/aimnetcentral
     - Elements: H, B, C, N, O, F, Si, P, S, Cl, As, Se, Br, I

- aimnet2nse
     - Description: Open-shell chemistry
     - Link: https://github.com/isayevlab/aimnetcentral
     - Elements: H, B, C, N, O, F, Si, P, S, Cl, As, Se, Br, I

- aimnet2pd
     - Description: Palladium-containing systems
     - Link: https://github.com/isayevlab/aimnetcentral
     - Elements: H, B, C, N, O, F, Si, P, S, Cl, Se, Br, Pd, I



- ani1ccx
     - Description: The ANI-1ccx model is an ensemble of 8 networks that was
       trained on the ANI-1ccx dataset, using transfer learning. The target
       accuracy is CCSD(T)*/CBS (CCSD(T) using the DPLNO-CCSD(T) method). It
       predicts energies on HCNO elements exclusively, it shouldn't be used
       with other atom types.
     - Link: https://github.com/aiqm/torchani
     - Link: https://doi.org/10.1038/s41597-020-0473-z
     - Elements: H, C, N, O
     
- ani1x
     - Description: The ANI-1x model is an ensemble of 8 networks that was
       trained using active learning on the ANI-1x dataset, the target level
       of theory is wB97X/6-31G(d). It predicts energies on HCNO elements
       exclusively, it shouldn't be used with other atom types.
     - Link: https://github.com/aiqm/torchani
     - Link: https://doi.org/10.1038/s41597-020-0473-z
     - Elements: H, C, N, O
     
- ani2x
     - Description: The ANI-2x model is an ensemble of 8 networks that was
       trained on the ANI-2x dataset. The target level of theory is
       wB97X/6-31G(d). It predicts energies on HCNOFSCl elements exclusively
       it shouldn't be used with other atom types.
     - Link: https://github.com/aiqm/torchani
     - Elements: H, C, N, O, F, S, Cl

- ani1xbb
     - Description: An ANI-Based Reactive Potential for Small Organic Molecules
     - Link: https://doi.org/10.1021/acs.jctc.5c00347
     - Elements: H, C, N, O
- pm6ml
     - Description: PM6-ML: The Synergy of Semiempirical Quantum Chemistry and
       Machine Learning Transformed into a Practical Computational Method
     - Link: https://doi.org/10.1021/acs.jctc.4c01330
     - Elements: H, C, N, O, P, S, F, Cl, Br, I, Li, Na, K, Mg, Ca

- fennix-bio1m
     - Description: A Foundation Model for Accurate Atomistic Simulations in Drug Design.
     - Description: Medium sized fit to DFT functional: wB97M-D3BJ / aug-ccpVTZ / ccECP
     - Link: https://doi.org/10.26434/chemrxiv-2025-f1hgn-v4
     - Link: https://doi.org/10.48550/arXiv.2405.01491
     - Link: https://github.com/FeNNol-tools/FeNNol-PMC
     - Elements: B, Br, C, Ca, Cl, F, H, I, K, Li, Mg, N, Na, O, P, S, Si, Zn

- fennix-bio1s
     - Description: A Foundation Model for Accurate Atomistic Simulations in Drug Design.
     - Description: Small sized fit to DFT functional: wB97M-D3BJ / aug-ccpVTZ / ccECP
     - Link: https://doi.org/10.26434/chemrxiv-2025-f1hgn-v4
     - Link: https://doi.org/10.48550/arXiv.2405.01491
     - Link: https://github.com/FeNNol-tools/FeNNol-PMC
     - Elements: B, Br, C, Ca, Cl, F, H, I, K, Li, Mg, N, Na, O, P, S, Si, Zn

- orb-v3-direct-inf-omat
     - Description: A "direct" model doesn't appear to calculate forces from backpropagation and may not conserve energy, but it is much faster than a conservative model.
     - Description: Trained to the OMol25 dataset (wB97M-V/def2-TZVPD)
     - Link: https://github.com/orbital-materials/orb-models
     - Link: https://arxiv.org/abs/2504.06231
     - Elements: Authors did not specify, but OMol25 contains 83 elements

- orb-v3-conservative-inf-omat
     - Description: Conservative models compute forces and stress via backpropagation, which is a physically motivated choice that appears necessary for certain types of simulation such as NVE Molecular dynamics. Conservative models are significantly slower and use more memory than their direct counterparts.
     - Description: Trained to the OMol25 dataset (wB97M-V/def2-TZVPD)
     - Link: https://github.com/orbital-materials/orb-models
     - Link: https://arxiv.org/abs/2504.06231
     - Elements: Authors did not specify, but OMol25 contains 83 elements
- am1
    - Description: AM1 method run through mopac
    - Link: 
    - Elements: H, Li, Be, B, C, N, O, F, Na, Mg, Al, Si, P, S, Cl, Cr, Zn, Ge, Br, Sn, I, Hg
    
- mndo
    - Description: AM1 method run through mopac
    - Link: 
    - Elements: H, Li, Be, B, C, N, O, F, Al, Si, P, S, Cl, Cr, Ge, Br, Sn, Hg, Pb, I
    
- mndod
    - Description: semiempirical method run through mopac
    - Link: 
    - Elements: H, Be, B, C, N, O, F, Na, Mg, Al, Si, P, S, Cl, Zn, Ge, Br, Cd, Sn, I, Hg, Pb
    
- pm3
    - Description: semiempirical method run through mopac
    - Link: 
    - Elements: H, Li, Be, B, C, N, O, F, Na, Mg, Al, Si, P, S, Cl, K, Ca, Sc, Ti, V, Cr, Mn, Fe, Co, Ni, Cu, Zn, Ga, Ge, As, Se, Br, Rb, Sr, Y, Zr, Nb, Mo, Tc, Ru, Rh, Pd, Ag, Cd, In, Sn, Sb, Te, I, Cs, Ba, Hg, Tl, Pb, Bi
    
- pm6
    - Description: semiempirical method run through mopac
    - Link: 
    - Elements: H, He, Li, Be, B, C, N, O, F, Ne, Na, Mg, Al, Si, P, S, Cl, Ar, K, Ca, Sc, Ti, V, Cr, Mn, Fe, Co, Ni, Cu, Zn, Ga, Ge, As, Se, Br, Kr, Rb, Sr, Y, Zr, Nb, Mo, Tc, Ru, Rh, Pd, Ag, Cd, In, Sn, Sb, Te, I, Xe, Cs, Ba, La, Lu, Hf, Ta, W, Re, Os, Ir, Pt, Au, Hg, Tl, Pb, Bi
    
- pm6-d3
    - Description: semiempirical method run through mopac
    - Link: 
    - Elements: H, He, Li, Be, B, C, N, O, F, Ne, Na, Mg, Al, Si, P, S, Cl, Ar, K, Ca, Sc, Ti, V, Cr, Mn, Fe, Co, Ni, Cu, Zn, Ga, Ge, As, Se, Br, Kr, Rb, Sr, Y, Zr, Nb, Mo, Tc, Ru, Rh, Pd, Ag, Cd, In, Sn, Sb, Te, I, Xe, Cs, Ba, La, Lu, Hf, Ta, W, Re, Os, Ir, Pt, Au, Hg, Tl, Pb, Bi, Po, At, Rn, Fr, Ra, Ac, Th, Pa, U, Np, Pu
    
- pm6-dh+
    - Description: semiempirical method run through mopac
    - Link: 
    - Elements: H, C, N, O, F, P, S, Cl, Br
    
- pm6-dh2
    - Description: semiempirical method run through mopac
    - Link: 
    - Elements: H, C, N, O, P, S, F, Cl, Br, I, Li, Na, K, Mg, Ca
    
- pm6-dh2x
    - Description: semiempirical method run through mopac
    - Link: 
    - Elements: H, He, Li, Be, B, C, N, O, F, Ne, Na, Mg, Al, Si, P, S, Cl, Ar, K, Ca, Sc, Ti, V, Cr, Mn, Fe, Co, Ni, Cu, Zn, Ga, Ge, As, Se, Br, Kr, Rb, Sr, Y, Zr, Nb, Mo, Tc, Ru, Rh, Pd, Ag, Cd, In, Sn, Sb, Te, I, Xe, Cs, Ba, La, Hf, Ta, W, Re, Os, Ir, Pt, Au, Hg, Tl, Pb, Bi, Lu
    
- pm6-d3h4
    - Description: semiempirical method run through mopac
    - Link: 
    - Elements: H, C, N, O, F, P, S, Cl, Br, I, Li, Na, K, Mg, Ca
    
- pm6-d3h4x
    - Description: semiempirical method run through mopac
    - Link: 
    - Elements: H, C, N, O, F, P, S, Cl, Br, I
    
- pmep
    - Description: semiempirical method run through mopac
    - Link: 
    - Elements: H, C, N, O, F, P, S, Cl, Br
    
- pm7
    - Description: semiempirical method run through mopac
    - Link: 
    - Elements: H, He, Li, Be, B, C, N, O, F, Ne, Na, Mg, Al, Si, P, S, Cl, Ar, K, Ca, Sc, Ti, V, Cr, Mn, Fe, Co, Ni, Cu, Zn, Ga, Ge, As, Se, Br, Kr, Rb, Sr, Y, Zr, Nb, Mo, Tc, Ru, Rh, Pd, Ag, Cd, In, Sn, Sb, Te, I, Xe, Cs, Ba, La, Lu, Hf, Ta, W, Re, Os, Ir, Pt, Au, Hg, Tl, Pb, Bi
    
- pm7-ts
    - Description: semiempirical method run through mopac
    - Link: 
    - Elements: H, He, Li, Be, B, C, N, O, F, Ne, Na, Mg, Al, Si, P, S, Cl, Ar, K, Ca, Sc, Ti, V, Cr, Mn, Fe, Co, Ni, Cu, Zn, Ga, Ge, As, Se, Br, Kr, Rb, Sr, Y, Zr, Nb, Mo, Tc, Ru, Rh, Pd, Ag, Cd, In, Sn, Sb, Te, I, Xe, Cs, Ba, La, Ce, Pr, Nd, Pm, Sm, Eu, Gd, Tb, Dy, Ho, Er, Tm, Yb, Lu, Hf, Ta, W, Re, Os, Ir, Pt, Au, Hg, Tl, Pb, Bi, Po, At, Rn, Fr, Ra, Ac, Th, Pa, U
    
- rm1
    - Description: semiempirical method run through mopac
    - Link: 
    - Elements: H, C, N, O, P, S, F, Cl, Br, and I
    

The following methods are available only if you:

1. Install a separate conda environment and install
   ffpopt with:
   ACADEMIC=TRUE python3 -m pip install --group fairchem .
2. Create an account on https://huggingface.co/
3. Visit https://huggingface.co/facebook/OMol25 and
   request permission from the maintainers to access
   the models via their online form. Then wait for
   access to be granted.
4. Visit https://huggingface.co/settings/tokens
   and create an access token. Write down a copy
   of the access token.
5. Use the "hf" command (installed within group fairchem)
   to authenticate by typing: "hf auth login"
   At the prompt, write your access token and press enter.
   Answer "n" when asked "Add token as git credential?"
6. You should see "Login successful".
7. You can then use the models listed below.

- OMOL25-ESEN-SM-DIRECT
     Description: 
     Link: https://huggingface.co/facebook/OMol25
     Link: https://arxiv.org/abs/2505.08762
     Elements: 

- OMOL25-ESEN-SM-CONSERVING
     Description: 
     Link: https://huggingface.co/facebook/OMol25
     Link: https://arxiv.org/abs/2505.08762
     Elements:
     
- OMOL25-ESEN-MD-DIRECT
     Description: 
     Link: https://huggingface.co/facebook/OMol25
     Link: https://arxiv.org/abs/2505.08762
     Elements:
     
- OMOL25-ESEN-LG-DIRECT
     Description: 
     Link: https://huggingface.co/facebook/OMol25
     Link: https://arxiv.org/abs/2505.08762
     Elements: 
    

The following methods are not available in the list of models,
but they are used by ffpopt-RespFit and ffpopt-DeltaRespFit
to preduct partial charges.

- espaloma
     Description: Conformation-independent ML network trained to reproduce AM1-BCC charges. Only available in pytorch version of ffpopt.
     Link: https://doi.org/10.1021/acs.jpca.4c01287
     Link: https://github.com/choderalab/espaloma-charge
     Elements: H,C,N,O,P,S,B,F,Cl,Br,I

- hilfiker
     Description: conformation-dependent ML network trained to reproduce PBE0-D4(BJ)/def2-TZVP. Only available in pytorch version of ffpopt.
     Link: https://doi.org/10.48550/arXiv.2512.13579
     Link: https://github.com/mathilfiker/ml_for_charges
     Elements: H,C,N,O,F,P,S,Cl
