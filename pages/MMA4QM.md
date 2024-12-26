---
layout: default
title: Mathematica for Quantum Mechanics
---

# Mathematica for Quantum Mechanics

Here we list examples of solving quantum mechanical problems with Mathematica.

## Change of basis for a spin-1/2 system

Below is the Q2.1 of the textbook: Schmied, R. (2020) Using Mathematica for Quantum Mechanics: A Student’s Manual. Singapore: Springer. Available at: https://doi.org/10.1007/978-981-13-7588-0.

**Q2.1** We describe a spin- $1 / 2$ system in the basis $\mathcal{B}$ containing the two states

$$
\begin{aligned}
& \left|\Uparrow_{\vartheta, \varphi}\right\rangle=\cos \left(\frac{\vartheta}{2}\right)|\uparrow\rangle+e^{\mathrm{i} \varphi} \sin \left(\frac{\vartheta}{2}\right)|\downarrow\rangle \\
& \left|\Downarrow_{\vartheta, \varphi}\right\rangle=-e^{-\mathrm{i} \varphi} \sin \left(\frac{\vartheta}{2}\right)|\uparrow\rangle+\cos \left(\frac{\vartheta}{2}\right)|\downarrow\rangle
\end{aligned}
$$

1. Show that the basis $\mathcal{B}=\left\{\left|\Uparrow_{\vartheta, \varphi}\right\rangle,\left|\Downarrow_{\vartheta, \varphi}\right\rangle\right\}$ is orthonormal.

2. Show that the basis $\mathcal{B}$ is complete: $\left|\Uparrow_{\vartheta, \varphi}\right\rangle\left\langle\Uparrow_{\vartheta, \varphi}\right|+\left|\Downarrow_{\vartheta, \varphi}\right\rangle\left\langle\Downarrow_{\vartheta, \varphi}\right|=\mathbb{1}$.

3. Express the states $|\uparrow\rangle$ and $|\downarrow\rangle$ as vectors in the basis $\mathcal{B}$.

4. Express the Pauli operators $\hat{\sigma}_x, \hat{\sigma}_y, \hat{\sigma}_z$ as matrices in the basis $\mathcal{B}$.

5. Show that $\left|\Uparrow_{\vartheta, \varphi}\right\rangle$ and $\left|\Downarrow_{\vartheta, \varphi}\right\rangle$ are eigenvectors of $\hat{\sigma}(\vartheta, \varphi)=\hat{\sigma}_x \sin (\vartheta) \cos (\varphi)+$ $\hat{\sigma}_y \sin (\vartheta) \sin (\varphi)+\hat{\sigma}_z \cos (\vartheta)$. What are the eigenvalues?

**Example solution**

```Mathematica
Q2.1

1.

In[149]:= up[\[Theta]_, \[CurlyPhi]_] := 
  Cos[\[Theta]/2]*{1, 0} + Exp[I*\[CurlyPhi]] Sin[\[Theta]/2]*{0, 1};
dn[\[Theta]_, \[CurlyPhi]_] := -Exp[-I*\[CurlyPhi]] Sin[\[Theta]/2]*{1, 0} + 
   Cos[\[Theta]/2]*{0, 1};

To prove that the basis {up, dn} is orthonormal, we calculate the following inner products:

In[151]:= Conjugate[up[\[Theta], \[CurlyPhi]]] . up[\[Theta], \[CurlyPhi]] // 
  ComplexExpand // FullSimplify
Conjugate[up[\[Theta], \[CurlyPhi]]] . dn[\[Theta], \[CurlyPhi]] // 
  ComplexExpand // FullSimplify
Conjugate[dn[\[Theta], \[CurlyPhi]]] . up[\[Theta], \[CurlyPhi]] // 
  ComplexExpand // FullSimplify
Conjugate[dn[\[Theta], \[CurlyPhi]]] . dn[\[Theta], \[CurlyPhi]] // 
  ComplexExpand // FullSimplify

Out[151]= 1

Out[152]= 0

Out[153]= 0

Out[154]= 1

2.

The proof of completeness goes straightforward:

In[155]:= KroneckerProduct[up[\[Theta], \[CurlyPhi]], 
    Conjugate[up[\[Theta], \[CurlyPhi]]]] + 
   KroneckerProduct[dn[\[Theta], \[CurlyPhi]], 
    Conjugate[dn[\[Theta], \[CurlyPhi]]]] // ComplexExpand // FullSimplify

Out[155]= {{1, 0}, {0, 1}}

3.

The coefficients are

In[156]:= Conjugate[up[\[Theta], \[CurlyPhi]]] . {1, 0} // 
  ComplexExpand // FullSimplify
Conjugate[dn[\[Theta], \[CurlyPhi]]] . {1, 0} // ComplexExpand // FullSimplify

Out[156]= Cos[\[Theta]/2]

Out[157]= -E^(I \[CurlyPhi]) Sin[\[Theta]/2]

4.

In[158]:= \[Sigma]x = PauliMatrix[1];
Conjugate[up[\[Theta], \[CurlyPhi]]] . \[Sigma]x . up[\[Theta], \[CurlyPhi]] //
   ComplexExpand // FullSimplify
Conjugate[up[\[Theta], \[CurlyPhi]]] . \[Sigma]x . dn[\[Theta], \[CurlyPhi]] //
   ComplexExpand // FullSimplify
Conjugate[dn[\[Theta], \[CurlyPhi]]] . \[Sigma]x . up[\[Theta], \[CurlyPhi]] //
   ComplexExpand // FullSimplify
Conjugate[dn[\[Theta], \[CurlyPhi]]] . \[Sigma]x . dn[\[Theta], \[CurlyPhi]] //
   ComplexExpand // FullSimplify

Out[159]= Cos[\[CurlyPhi]] Sin[\[Theta]]

Out[160]= E^(-I \[CurlyPhi]) (Cos[\[Theta]] Cos[\[CurlyPhi]] + I Sin[\[CurlyPhi]])

Out[161]= E^(I \[CurlyPhi]) (Cos[\[Theta]] Cos[\[CurlyPhi]] - I Sin[\[CurlyPhi]])

Out[162]= -Cos[\[CurlyPhi]] Sin[\[Theta]]

Processes for other Pauli matrices go similar.

5.

In[163]:= \[Sigma]y = PauliMatrix[2];
\[Sigma]z = PauliMatrix[3];
\[Sigma][\[Theta]_, \[CurlyPhi]_] := 
 Sin[\[Theta]] Cos[\[CurlyPhi]] \[Sigma]x + 
  Sin[\[Theta]] Sin[\[CurlyPhi]] \[Sigma]y + Cos[\[Theta]] \[Sigma]z
Eigenvalues[\[Sigma][\[Theta], \[CurlyPhi]]]
\[Sigma][\[Theta], \[CurlyPhi]] . up[\[Theta], \[CurlyPhi]] == 
  up[\[Theta], \[CurlyPhi]] // FullSimplify
\[Sigma][\[Theta], \[CurlyPhi]] . 
   dn[\[Theta], \[CurlyPhi]] == -dn[\[Theta], \[CurlyPhi]] // FullSimplify

Out[166]= {-1, 1}

Out[167]= True

Out[168]= True
```