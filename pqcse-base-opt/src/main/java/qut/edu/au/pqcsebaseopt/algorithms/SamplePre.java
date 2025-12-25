package qut.edu.au.pqcsebaseopt.algorithms;

import qut.edu.au.pqcsebaseopt.tools.BigIntMatrix;

import java.math.BigInteger;
import java.util.Arrays;
import static qut.edu.au.pqcsebaseopt.algorithms.SampleD.sampleD;

public class SamplePre {

    public static BigInteger[] samplePre(BigIntMatrix A, BigIntMatrix B, double sigma, BigInteger[] u) {
        int m = A.getCols();

        // 1. Solve A*t = u using Optimized Solver
        BigInteger[] t = solveLinearSystemFast(A, u);

        if (t == null) {
            throw new RuntimeException("Linear system A*t = u has no solution.");
        }

        // 2. Prepare center c = -t
        double[] c = new double[m];
        for (int i = 0; i < m; i++) {
            c[i] = -t[i].doubleValue();
        }

        // 3. Sample D.
        // Note: For maximum performance, B'sigma GSO should be precomputed and passed in.
        // Here we use the optimized internal logic of SampleD.
        BigInteger[] x = sampleD(B, sigma, c, A.getRows());

        // 4. v = x + t
        BigInteger[] v = new BigInteger[m];
        Arrays.setAll(v, i -> x[i].add(t[i]));

        return v;
    }

    /**
     * [Optimization] Fast SamplePre with Precomputed GSO.
     * Prevents OutOfMemoryError by reusing GSO data across parallel threads.
     */
    public static BigInteger[] samplePreFast(BigIntMatrix A, BigIntMatrix B_basis,
                                             double[][] B_tilde, double[] normsSq, double[][] B_double,
                                             double sigma, BigInteger[] u) {
        int m = A.getCols();

        // 1. Solve Linear System
        BigInteger[] t = solveLinearSystemFast(A, u);
        if (t == null) throw new RuntimeException("No solution");

        // 2. Center c = -t
        double[] c = new double[m];
        for (int i = 0; i < m; i++) c[i] = -t[i].doubleValue();

        // 3. Sample D using fast GSO (No re-allocation!)
        BigInteger[] x = qut.edu.au.pqcsebaseopt.algorithms.SampleD.sampleD_Fast(
                B_tilde, normsSq, B_double, B_basis, sigma, c, A.getRows());

        // 4. v = x + t
        BigInteger[] v = new BigInteger[m];
        Arrays.setAll(v, i -> x[i].add(t[i]));
        return v;
    }

    /**
     * Optimized Linear System Solver.
     * Strategy: Try to invert the first n x n block first.
     * Since A = [A1 | A2] and A1 is usually random/invertible, this works 99% of the time.
     * Cost: O(n^3) vs O(n^2 * m) for full RREF.
     */
    public static BigInteger[] solveLinearSystemFast(BigIntMatrix A, BigInteger[] u) {
        int n = A.getRows();
        int m = A.getCols();
        BigInteger q = A.getModulus();

        // 1. Try to extract A1 (first n columns)
        try {
            BigIntMatrix A1 = A.getSubMatrix(0, n - 1, 0, n - 1);
            // Compute A1^-1 mod q
            BigIntMatrix A1_inv = BigIntMatrix.inverseQ(A1); // O(n^3)

            // If successful, t_head = A1^-1 * u
            // t_tail = 0
            BigInteger[] t_head = BigIntMatrix.multiplyMatrixVector(A1_inv, u);

            BigInteger[] t = new BigInteger[m];
            System.arraycopy(t_head, 0, t, 0, n);
            for (int i = n; i < m; i++) t[i] = BigInteger.ZERO;

            return t;

        } catch (RuntimeException e) {
            // Fallback: A1 was singular, run full Gaussian elimination
            // System.out.println("Fast solve failed, falling back to RREF...");
            return solveLinearSystemRREF(A, u);
        }
    }

    /**
     * Fallback RREF Solver (Original Logic)
     */
    public static BigInteger[] solveLinearSystemRREF(BigIntMatrix A, BigInteger[] b) {
        // ... (保持你原有的 Gaussian elimination 代码不变作为 fallback) ...
        // 为了节省篇幅，这里复用你原有的代码逻辑
        // 请保留原本的 solveLinearSystem 实现并重命名为 solveLinearSystemRREF
        int n = A.getRows();
        int m = A.getCols();
        BigInteger q = A.getModulus();

        BigInteger[][] M = new BigInteger[n][m + 1];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                M[i][j] = A.get(i, j).mod(q);
            }
            M[i][m] = b[i].mod(q);
        }

        int pivotRow = 0;
        int[] pivotColIndex = new int[n];
        Arrays.fill(pivotColIndex, -1);

        for (int col = 0; col < m && pivotRow < n; col++) {
            int sel = -1;
            for (int row = pivotRow; row < n; row++) {
                if (!M[row][col].gcd(q).equals(BigInteger.ONE)) {
                    if (!q.isProbablePrime(20) && !M[row][col].equals(BigInteger.ZERO)) continue;
                    else if (q.isProbablePrime(20) && !M[row][col].equals(BigInteger.ZERO)) { sel = row; break; }
                } else {
                    sel = row; break;
                }
            }
            if (sel == -1) continue;

            if (sel != pivotRow) {
                BigInteger[] temp = M[pivotRow]; M[pivotRow] = M[sel]; M[sel] = temp;
            }
            pivotColIndex[pivotRow] = col;
            BigInteger pivot = M[pivotRow][col];
            BigInteger inv = pivot.modInverse(q);
            for (int j = col; j <= m; j++) M[pivotRow][j] = M[pivotRow][j].multiply(inv).mod(q);
            for (int i = 0; i < n; i++) {
                if (i != pivotRow) {
                    BigInteger factor = M[i][col];
                    if (!factor.equals(BigInteger.ZERO)) {
                        for (int j = col; j <= m; j++) M[i][j] = M[i][j].subtract(factor.multiply(M[pivotRow][j])).mod(q);
                    }
                }
            }
            pivotRow++;
        }

        BigInteger[] x = new BigInteger[m];
        Arrays.fill(x, BigInteger.ZERO);
        for (int i = n - 1; i >= 0; i--) {
            int pCol = pivotColIndex[i];
            if (pCol != -1) x[pCol] = M[i][m];
        }
        return x;
    }

    /**
     * Calculates the maximum length of the Gram-Schmidt orthogonalized vectors (||S~||).
     * This value is crucial for determining the lower bound of the Gaussian parameter sigma.
     *
     * @param B The lattice basis matrix. **Assumes Column Basis** (each column is a vector).
     * @return The maximum Euclidean norm of the orthogonalized basis vectors.
     */
    public static double computeGaussParam(BigIntMatrix B) {
        int n = B.getCols(); // Number of basis vectors
        int m = B.getRows();    // Dimension of vectors

        // 2. Extract data to double array for GSO computation
        // B_double[i] represents the i-th basis vector (which is the i-th COLUMN of B)
        double[][] B_double = new double[n][m];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                // Read the i-th column, j-th row
                B_double[i][j] = B.get(j, i).doubleValue();
            }
        }

        // 3. Perform Gram-Schmidt Orthogonalization
        double[][] B_tilde = new double[n][m];
        double[] normsSq = new double[n]; // Stores squared norms ||b_tilde_i||^2

        for (int i = 0; i < n; i++) {
            // Initialize: B_tilde[i] = B[i]
            System.arraycopy(B_double[i], 0, B_tilde[i], 0, m);

            // Subtract projections onto all previous orthogonal vectors
            for (int j = 0; j < i; j++) {
                // Compute dot product <B[i], B_tilde[j]>
                double dotProd = 0;
                for (int k = 0; k < m; k++) {
                    dotProd += B_double[i][k] * B_tilde[j][k];
                }

                // coefficient = <B[i], B_tilde[j]> / ||B_tilde[j]||^2
                double coef = dotProd / normsSq[j];

                // B_tilde[i] = B_tilde[i] - coef * B_tilde[j]
                for (int k = 0; k < m; k++) {
                    B_tilde[i][k] -= coef * B_tilde[j][k];
                }
            }

            // Compute and store the squared norm of the orthogonalized vector
            double norm = 0;
            for (int k = 0; k < m; k++) {
                norm += B_tilde[i][k] * B_tilde[i][k];
            }
            normsSq[i] = norm;
        }

        // 4. Find the maximum squared norm
        double maxNormSq = 0;
        for (double val : normsSq) {
            if (val > maxNormSq) {
                maxNormSq = val;
            }
        }

        // 5. Return the actual maximum length (Square Root)
        return Math.sqrt(maxNormSq);
    }
}