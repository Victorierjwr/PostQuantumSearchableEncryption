package qut.edu.au.pqcsebaseopt.tools;

import java.math.BigInteger;
import java.util.stream.IntStream;

/**
 * Optimized Hermite Normal Form calculation.
 * Adapts to the high-performance 1D BigIntMatrix.
 * Uses parallel streams for vector operations to handle N=4093.
 */
public class HermiteNormalForm {

    /**
     * Compute HNF for the optimized BigIntMatrix.
     * Note: HNF intermediate values GROW exponentially.
     * We convert primitive longs to BigIntegers internally for calculation.
     */
    public static BigIntMatrix hermiteNormalForm(BigIntMatrix input) {
        int rows = input.getRows();
        int cols = input.getCols();

        // 1. 转置并提升为 BigInteger (HNF core logic works on transpose)
        // Flattened 1D array for Transpose: T_rows = cols, T_cols = rows
        BigInteger[] workMatrix = new BigInteger[cols * rows];

        // Parallel Transpose & Expansion (long -> BigInteger)
        // input(i, j) -> work(j, i)
        IntStream.range(0, rows).parallel().forEach(i -> {
            for (int j = 0; j < cols; j++) {
                // getLong avoids BigInteger creation during fetch
                long val = input.getLong(i, j);
                workMatrix[j * rows + i] = BigInteger.valueOf(val);
            }
        });

        // 2. Compute Lower Triangular HNF on the transposed matrix
        // (Working dimensions are now m=cols, n=rows)
        hnfLowerParallel(workMatrix, cols, rows);

        // 3. Transpose back and pack into result BigIntMatrix
        // Result dimensions: rows x cols
        BigIntMatrix res = new BigIntMatrix(rows, cols, input.getModulus());

        // Since HNF result might exceed long (mod q), usually we expect
        // the result in TrapGen to be used Mod q eventually.
        // But strictly, HNF returns integers.
        // If result values > 2^63, BigIntMatrix will truncate if we use set(long).
        // So we must use set(BigInteger).

        IntStream.range(0, cols).parallel().forEach(i -> {
            for (int j = 0; j < rows; j++) {
                // work(i, j) -> res(j, i)
                BigInteger val = workMatrix[i * rows + j];
                // Warning: If val > modulus, this preserves the integer value
                // but optimized BigIntMatrix might treat it as long internally if < 2^63.
                res.set(j, i, val);
            }
        });

        return res;
    }

    /**
     * Parallelized Lower HNF Algorithm.
     * Operates on a flattened 1D array 'M' of size m*n.
     */
    private static void hnfLowerParallel(BigInteger[] M, int m, int n) {
        // Iterate columns from right to left
        for (int c = n - 1; c >= 0; c--) {
            int pivotRow = -1;

            // Find pivot
            int diagIdx = c * n + c; // M[c][c]
            if (c < m && !M[diagIdx].equals(BigInteger.ZERO)) {
                pivotRow = c;
            } else {
                for (int i = Math.min(c, m - 1); i >= 0; i--) {
                    if (!M[i * n + c].equals(BigInteger.ZERO)) {
                        pivotRow = i;
                        break;
                    }
                }
            }

            if (pivotRow == -1) continue;

            // Swap rows to place pivot
            if (pivotRow != c && c < m) {
                swapRows(M, c, pivotRow, n);
            }
            int currPivotRow = (c < m) ? c : pivotRow;

            // Eliminate above pivot (Rows 0 to c-1)
            // This loop is sequential (iterative elimination),
            // but the row operations inside are PARALLEL.
            for (int i = c - 1; i >= 0; i--) {
                int idxPivot = currPivotRow * n + c;
                int idxTarget = i * n + c;

                if (M[idxTarget].equals(BigInteger.ZERO)) continue;

                BigInteger u = M[idxPivot];
                BigInteger v = M[idxTarget];

                // Standard Extended GCD
                BigInteger[] eg = extendedGcdBig(u, v);
                BigInteger g = eg[0];
                BigInteger x = eg[1];
                BigInteger y = eg[2];

                BigInteger uDivG = u.divide(g);
                BigInteger vDivG = v.divide(g);
                BigInteger negVDivG = vDivG.negate();

                // Apply Row Operations in PARALLEL
                // Row_Pivot = x * Row_Pivot + y * Row_Target
                // Row_Target = -v/g * Row_Pivot + u/g * Row_Target
                // For N=4093, this parallelism is crucial.
                final int pRow = currPivotRow;
                final int tRow = i;

                IntStream.range(0, n).parallel().forEach(k -> {
                    int pIdx = pRow * n + k;
                    int tIdx = tRow * n + k;

                    BigInteger valP = M[pIdx];
                    BigInteger valI = M[tIdx];

                    // Calculation involves creating new BigIntegers, unavoidable in HNF
                    M[pIdx] = valP.multiply(x).add(valI.multiply(y));
                    M[tIdx] = valP.multiply(negVDivG).add(valI.multiply(uDivG));
                });
            }

            // Ensure pivot is positive
            int pivotIdxFinal = currPivotRow * n + c;
            if (M[pivotIdxFinal].signum() < 0) {
                final int rowToFlip = currPivotRow;
                IntStream.range(0, n).parallel().forEach(k -> {
                    int idx = rowToFlip * n + k;
                    M[idx] = M[idx].negate();
                });
            }

            // Reduce below pivot (Rows c+1 to m-1)
            for (int i = c + 1; i < m; i++) {
                int idxT = i * n + c;
                if (M[idxT].equals(BigInteger.ZERO)) continue;

                BigInteger pivot = M[currPivotRow * n + c];
                BigInteger val = M[idxT];

                BigInteger q = val.divide(pivot);
                // Floor adjustment for negative division
                if (val.signum() < 0 && !val.mod(pivot).equals(BigInteger.ZERO)) {
                    q = q.subtract(BigInteger.ONE);
                }

                if (!q.equals(BigInteger.ZERO)) {
                    BigInteger negQ = q.negate();
                    final int targetRow = i;
                    final int srcRow = currPivotRow;

                    // Parallel Row Combination: Target += src * (-q)
                    IntStream.range(0, n).parallel().forEach(k -> {
                        int sIdx = srcRow * n + k;
                        int tIdx = targetRow * n + k;
                        M[tIdx] = M[tIdx].add(negQ.multiply(M[sIdx]));
                    });
                }
            }
        }
    }

    // --- Utilities ---

    private static void swapRows(BigInteger[] M, int r1, int r2, int cols) {
        // Parallel Swap
        IntStream.range(0, cols).parallel().forEach(k -> {
            int idx1 = r1 * cols + k;
            int idx2 = r2 * cols + k;
            BigInteger tmp = M[idx1];
            M[idx1] = M[idx2];
            M[idx2] = tmp;
        });
    }

    private static BigInteger[] extendedGcdBig(BigInteger a, BigInteger b) {
        BigInteger x = BigInteger.ZERO, y = BigInteger.ONE;
        BigInteger lastX = BigInteger.ONE, lastY = BigInteger.ZERO;
        BigInteger temp;
        while (!b.equals(BigInteger.ZERO)) {
            BigInteger[] dr = a.divideAndRemainder(b);
            BigInteger q = dr[0];
            BigInteger r = dr[1];
            a = b; b = r;
            temp = x; x = lastX.subtract(q.multiply(x)); lastX = temp;
            temp = y; y = lastY.subtract(q.multiply(y)); lastY = temp;
        }
        if (a.signum() < 0) {
            a = a.negate(); lastX = lastX.negate(); lastY = lastY.negate();
        }
        return new BigInteger[]{a, lastX, lastY};
    }
}