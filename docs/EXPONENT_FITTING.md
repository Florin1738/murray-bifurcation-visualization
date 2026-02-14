# Murray Exponent Fitting - Cost Function Selection

## Date: 2026-01-21

## Problem Statement

The original implementation was finding an optimal Murray exponent of α=1.5 for all samples, which was highly suspicious and indicated a fundamental issue with the optimization approach.

## Root Cause

The original cost function was:

```python
minimize Σ(r_p^α - (r_c1^α + r_c2^α))²
```

**Critical Flaw:** This is **not scale-invariant**.

- Large vessels contribute exponentially more to the sum than small vessels
- A bifurcation with r_p=100 contributes ~10,000× more than one with r_p=10 (when α≈3)
- The optimizer minimized errors for the largest bifurcations only
- The fitted exponent did not reflect the true scaling relationship across the entire vascular network

## Solution: Scale-Invariant Relative Residuals

The new cost function normalizes each residual by the parent vessel's capacity:

```python
minimize Σ[(r_p^α - (r_c1^α + r_c2^α))² / (r_p^α)²]
```

Which simplifies to:

```python
minimize Σ[1 - phi(α)]²
```

where `phi(α) = (r_c1^α + r_c2^α) / r_p^α`

### Why This Works

1. **Scale Invariance**: Each bifurcation contributes equally regardless of vessel size
2. **Physical Meaning**: Minimizes relative deviations from the Murray branching rule
3. **Statistical Validity**: Treats all bifurcations as equally informative about the scaling relationship
4. **Ratio-Based**: Directly optimizes the phi ratio to be close to 1.0 across all scales

## Alternative Cost Functions Considered

### 1. Log-Space Residuals
```python
minimize Σ[log(r_p^α) - log(r_c1^α + r_c2^α)]²
```
- Also scale-invariant (multiplicative errors)
- More sensitive to small vessels
- Less robust to measurement noise at small scales
- **Not chosen** because we want equal weighting, not small-vessel emphasis

### 2. Absolute Residuals (original)
```python
minimize Σ(r_p^α - (r_c1^α + r_c2^α))²
```
- **Not scale-invariant**
- Dominated by large vessels
- **Rejected** - this was the source of the bug

### 3. Weighted by Vessel Importance
```python
minimize Σ w_i × [residual_i]²
```
where `w_i` could be flow, volume, or other importance metric
- Adds complexity and assumptions
- Requires additional flow modeling
- **Not chosen** - we want to characterize the geometric scaling first

## Theoretical Background

### Murray's Law and the Exponent α

In the general form, Murray's law states:

```
r_p^α = r_c1^α + r_c2^α
```

The exponent α encodes how "transport capacity" scales with radius in an optimized network:

- **α = 3**: Classical Murray's law for Newtonian flow with volume minimization
- **α = 2.5**: Surface area minimization constraint (Revellin et al.)
- **α = 2-3**: Range observed for blood vessels accounting for non-Newtonian effects
- **α > 3**: Structural/developmental constraints favor larger parent vessels
- **α < 3**: Metabolic/flow constraints favor more distributed capacity

### Flow-Radius Relationship

The branching exponent is related to the flow-radius scaling:

```
Q ∝ r^α
```

At a bifurcation with flow conservation:
```
Q_p = Q_c1 + Q_c2
k·r_p^α = k·r_c1^α + k·r_c2^α
r_p^α = r_c1^α + r_c2^α
```

The constant k cancels, yielding the branching rule.

## Implementation Details

The scale-invariant cost function is implemented in `fit_murray_exponent()` at line ~1061 in `skeleton_analysis_plots.py`.

### Search Range
- Bounds: [1.0, 5.0]
- Wide enough to capture biological variation
- Includes classical Murray (3.0) and surface-area optimized (2.5)

### Diagnostic Output
The function now prints:
- Bifurcation count and radii distributions
- Traditional Murray phi (α=3) statistics
- Objective function values at test points (α = 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0)
- Optimization success status

## Expected Results

With the corrected cost function, we expect:

1. **Different exponents across samples**: Each species/sample should have distinct optimal α
2. **Biologically plausible range**: α ∈ [2.0, 3.5] for most vascular networks
3. **Variation with biological context**: Different organs, species, or developmental stages may show different optima
4. **Not all α=1.5**: The bug that caused uniform α=1.5 should be resolved

## References

- Murray, C. D. (1926). The physiological principle of minimum work.
- Sherman, T. F. (1981). On connecting large vessels to small.
- Revellin, R., et al. (2009). Extension of Murray's law using a non-Newtonian model of blood flow.

## Verification

To verify the fix is working:
1. Run the script and check the diagnostic output
2. Confirm different samples show different optimal α values
3. Check that objective function values show a clear minimum (not monotonic)
4. Verify α values fall in biologically plausible range [2.0-3.5]
