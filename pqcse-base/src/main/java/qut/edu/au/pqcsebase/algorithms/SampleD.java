package qut.edu.au.pqcsebase.algorithms;

import qut.edu.au.pqcsebase.tools.BigIntMatrix;

import java.math.BigInteger;
import java.util.Arrays;

import static qut.edu.au.pqcsebase.algorithms.SampleZ.sampleZ;

/**
 * Gaussian distribution
 * Implementation of Discrete Gaussian Sampling over Lattices (SampleD / Randomized Nearest Plane).
 * References: GPV08, AP11 papers.
 */

public class SampleD {

    /**
     * Sample from the discrete Gaussian distribution over the lattice defined by the basis B,
     * with parameter sigma and center c.
     * also called SampleGaussian(A, S, \sigma, 0)$
     *
     * @param B The basis matrix of the lattice.Each row is a basis vector
     * @param sigma The Gaussian parameter.
     * @param c The center vector.
     * @param nParam The dimension of the lattice.
     * @return v sample from the discrete Gaussian distribution over the lattice.
     */
    public static BigInteger[] sampleD(BigIntMatrix B, double sigma, double[] c, long nParam) {
        int n = B.getColumnDimension();
        int m = B.getRowDimension();

        //for GSO computation
        double[][] B_double = new double[n][m];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                B_double[i][j] = B.get(j, i).doubleValue();
            }
        }

        //Gram-Schmidt
        double[][] B_tilde = new double[n][m];
        //B_tilde[i]||^2
        double[] normsSq = new double[n];

        //need long time
        computeGramSchmidt(B_double, B_tilde, normsSq);

        // Initialize v as zero vector
        BigInteger[] v = new BigInteger[m];
        Arrays.fill(v, BigInteger.ZERO);

        // Sampling process
        double[] c_current = Arrays.copyOf(c, c.length);
        for (int i = n - 1; i >= 0; i--) {
            // Compute mu_i = <c_current, B_tilde[i]> / ||B_tilde[i]||^2
            double c_prime = doProduct(c_current, B_tilde[i]) / normsSq[i];
            double s_prime = sigma / Math.sqrt(normsSq[i]);

            BigInteger z_i = sampleZ(s_prime, c_prime, nParam);
            double z_double = z_i.doubleValue();

            // c(i-1) = c(i) - z_i * B[i]
            for (int k = 0; k < m; k++) {
                c_current[k] = c_current[k] - z_double * B_double[i][k];
            }

            // v(i-1) = v(i) + z_i * B[i]
            for (int k = 0; k < m; k++) {
                BigInteger val = B.get(k, i).multiply(z_i);
                v[k] = v[k].add(val);
            }
        }
        return v;
    }

    /**
     * Helper: Computes Gram-Schmidt Orthogonalization.
     *
     * @param B       Input basis (rows are vectors).
     * @param B_tilde Output orthogonalized basis.
     * @param normsSq Output squared norms of orthogonal vectors.
     */
    public static void computeGramSchmidt(double[][] B, double[][] B_tilde, double[] normsSq) {
        int n = B.length;
        int m = B[0].length;

        for (int i = 0; i < n; i++) {
            // Initialize B_tilde[i] = B[i]
            System.arraycopy(B[i], 0, B_tilde[i], 0, m);

            // Subtract projections onto all previous orthogonal vectors
            for (int j = 0; j < i; j++) {
                // coef = <B[i], B_tilde[j]> / ||B_tilde[j]||^2
                double coef = doProduct(B[i], B_tilde[j]) / normsSq[j];

                // B_tilde[i] = B_tilde[i] - coef * B_tilde[j]
                for (int k = 0; k < m; k++) {
                    B_tilde[i][k] -= coef * B_tilde[j][k];
                }
            }

            // Compute and store squared norm
            normsSq[i] = doProduct(B_tilde[i], B_tilde[i]);
        }
    }

    /**
     * Helper: Computes dot product of two vectors.
     */
    private static double doProduct(double[] a, double[] b) {
        double sum = 0;
        for (int i = 0; i < a.length; i++) {
            sum += a[i] * b[i];
        }
        return sum;
    }

    /**
     * Fast version of SampleD when GSO is precomputed.
     *
     * @param B_tilde     Precomputed orthogonalized basis.
     * @param normsSq     Precomputed squared norms of orthogonal vectors.
     * @param B_original  Original basis matrix.
     * @param sigma           Gaussian parameter.
     * @param c           Center vector.
     * @param nParam      Security parameter.
     * @return v Sample from the discrete Gaussian distribution over the lattice.
     */
    public static BigInteger[] sampleD_Fast(double[][] B_tilde, double[] normsSq, BigIntMatrix B_original, double sigma, double[] c, long nParam) {
        int m = B_original.getRowDimension();
        int n = B_original.getColumnDimension();

        // 1. Initialization Phase
        BigInteger[] v = new BigInteger[m];
        Arrays.fill(v, BigInteger.ZERO);

        double[] c_current = Arrays.copyOf(c, c.length);

        for (int i = n - 1; i >= 0; i--) {
            double c_prime = doProduct(c_current, B_tilde[i]) / normsSq[i];
            double s_prime = sigma / Math.sqrt(normsSq[i]);

            BigInteger z_i = sampleZ(s_prime, c_prime, nParam);
            double z_double = z_i.doubleValue();

            for (int k = 0; k < m; k++) {
                double original_val = B_original.get(k, i).doubleValue();
                c_current[k] = c_current[k] - z_double * original_val;
            }

            for (int k = 0; k < m; k++) {
                BigInteger val = B_original.get(k, i).multiply(z_i);
                v[k] = v[k].add(val);
            }
        }
        return v;
    }
}
