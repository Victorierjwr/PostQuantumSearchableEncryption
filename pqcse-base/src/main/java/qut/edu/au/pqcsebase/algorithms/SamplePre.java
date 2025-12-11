package qut.edu.au.pqcsebase.algorithms;

import qut.edu.au.pqcsebase.tools.BigIntMatrix;

import java.math.BigInteger;
import java.util.Arrays;

import static qut.edu.au.pqcsebase.algorithms.SampleD.sampleD;

/**
 * preimage sampleable functions (PSFs).
 */
public class SamplePre {

    /**
     * SamplePre: Finds a short vector v such that A*v = u mod q.
     *
     * @param A The matrix A (n x m).
     * @param B The trapdoor basis matrix for the lattice Lambda_perp(A). TrapGen need to inverse B = TA^T
     * @param s The Gaussian parameter sigma.
     * @param u The target vector (length n).
     * @return A short vector v (length m) satisfying the equation.
     */
    public static BigInteger[] samplePre(BigIntMatrix A, BigIntMatrix B, double s, BigInteger[] u) {
        int n = A.getRowDimension();
        int m = A.getColumnDimension();
        BigInteger q = A.getModulus();

        // 1. Solve the linear system A * t = u mod q to find a particular solution t.
        BigInteger[] t = solveLinearSystem(A, u);

        if (t == null) {
            throw new RuntimeException("Linear system A*t = u has no solution (A might not be full rank).");
        }

        // 2. center c = -t for SampleD.
        double[] c = new double[m];
        for (int i = 0; i < m; i++) {
            // c[i] = -t[i] (converted to double)
            c[i] = -t[i].doubleValue();
        }

        // 3. discrete Gaussian sampling on the lattice: x ~ D_{Lambda_perp, s, -t}.
        BigInteger[] x = sampleD(B, s, c, n);

        // 4. Compute the final result: v = x + t.
        BigInteger[] v = new BigInteger[m];
        for (int i = 0; i < m; i++) {
            v[i] = x[i].add(t[i]);
        }

        return v;
    }

    /**
     * Solves the linear system A * x = b mod q using Gaussian elimination (Gauss-Jordan).
     * Since A is n x m (with m > n) and typically full rank, the system is under-determined.
     * We return one particular solution (usually by setting free variables to 0).
     *
     * @param A The coefficient matrix.
     * @param b The target vector.
     * @return A particular solution vector x, or null if no solution exists.
     */
    public static BigInteger[] solveLinearSystem(BigIntMatrix A, BigInteger[] b) {
        int n = A.getRowDimension();
        int m = A.getColumnDimension();
        BigInteger q = A.getModulus();

        // Construct the augmented matrix M = [A | b]
        // Dimension: n x (m + 1)
        BigInteger[][] M = new BigInteger[n][m + 1];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                M[i][j] = A.get(i, j).mod(q);
            }
            M[i][m] = b[i].mod(q);
        }

        // Gaussian elimination to Reduced Row Echelon Form (RREF)
        int pivotRow = 0;
        int[] pivotColIndex = new int[n]; // Tracks the column index of the pivot for each row
        Arrays.fill(pivotColIndex, -1);

        for (int col = 0; col < m && pivotRow < n; col++) {
            // Find a non-zero pivot in the current column
            int sel = -1;
            for (int row = pivotRow; row < n; row++) {
                // We prefer a pivot coprime to q to ensure invertibility (modular inverse exists).
                // If q is prime, any non-zero element works.
                if (!M[row][col].gcd(q).equals(BigInteger.ONE)) {
                    // If q is not prime, we skip rows where gcd(val, q) != 1 unless unavoidable.
                    // For prime q (typical in LWE), this check simplifies to !val.equals(ZERO).
                    if (!q.isProbablePrime(20) && !M[row][col].equals(BigInteger.ZERO)) {
                        continue;
                    } else if (q.isProbablePrime(20) && !M[row][col].equals(BigInteger.ZERO)) {
                        sel = row;
                        break;
                    }
                } else {
                    sel = row;
                    break;
                }
            }

            if (sel == -1) continue; // Column is all zeros (or non-invertible), skip to next column

            // Swap rows
            if (sel != pivotRow) {
                BigInteger[] temp = M[pivotRow];
                M[pivotRow] = M[sel];
                M[sel] = temp;
            }

            // Record pivot position
            pivotColIndex[pivotRow] = col;

            // Normalize the pivot row (make pivot = 1)
            BigInteger pivot = M[pivotRow][col];
            BigInteger inv = pivot.modInverse(q);
            for (int j = col; j <= m; j++) {
                M[pivotRow][j] = M[pivotRow][j].multiply(inv).mod(q);
            }

            // Eliminate other rows
            for (int i = 0; i < n; i++) {
                if (i != pivotRow) {
                    BigInteger factor = M[i][col];
                    if (!factor.equals(BigInteger.ZERO)) {
                        for (int j = col; j <= m; j++) {
                            M[i][j] = M[i][j].subtract(factor.multiply(M[pivotRow][j])).mod(q);
                        }
                    }
                }
            }
            pivotRow++;
        }

        // Back Substitution (Extraction of solution)
        // Initialize solution vector x to all zeros (treating free variables as 0)
        BigInteger[] x = new BigInteger[m];
        Arrays.fill(x, BigInteger.ZERO);

        for (int i = n - 1; i >= 0; i--) {
            int pCol = pivotColIndex[i];
            if (pCol != -1) {
                BigInteger val = M[i][m]; // The value in the augmented column
                // Since the matrix is in RREF (diagonalized at pivots), and free vars are 0,
                // x[pivot_col] is simply the value in the augmented column.
                x[pCol] = val;
            }
        }

        return x;
    }
}
