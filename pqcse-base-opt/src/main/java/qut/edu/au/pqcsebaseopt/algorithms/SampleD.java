package qut.edu.au.pqcsebaseopt.algorithms;


import qut.edu.au.pqcsebaseopt.tools.BigIntMatrix;

import java.math.BigInteger;
import java.util.Arrays;
import java.util.stream.IntStream;

import static qut.edu.au.pqcsebaseopt.algorithms.SampleZ.sampleZ;

public class SampleD {

    /**
     * Standard SampleD.
     * WARNING: This computes GSO on every call, which is very slow.
     * Use sampleD_Fast with precomputed GSO for better performance.
     */
    public static BigInteger[] sampleD(BigIntMatrix B, double sigma, double[] c, long nParam) {
        int n = B.getCols();
        int m = B.getRows();

        // 1. Convert B to double once
        double[][] B_double = new double[n][m];
        // Parallel conversion
        IntStream.range(0, n).parallel().forEach(i -> {
            for (int j = 0; j < m; j++) {
                B_double[i][j] = B.get(j, i).doubleValue();
            }
        });

        // 2. Compute GSO (Expensive!)
        double[][] B_tilde = new double[n][m];
        double[] normsSq = new double[n];
        computeGramSchmidt(B_double, B_tilde, normsSq);

        // 3. Delegate to fast implementation
        return sampleD_Fast(B_tilde, normsSq, B_double, B, sigma, c, nParam);
    }

    /**
     * Optimized GSO Computation (Parallelized outer loops where possible)
     */
    public static void computeGramSchmidt(double[][] B, double[][] B_tilde, double[] normsSq) {
        int n = B.length;
        int m = B[0].length;
        System.out.println("Start GSO (Optimized Sequential)");
        long startTime = System.currentTimeMillis();

        for (int i = 0; i < n; i++) {
            // 打印进度
            if (i % 100 == 0) {
                long currentTime = System.currentTimeMillis();
                System.out.println("GSO Iteration " + i + " of " + n + " (" + (currentTime - startTime)/1000 + "s)");
            }

            // 初始化: B_tilde[i] = B[i]
            System.arraycopy(B[i], 0, B_tilde[i], 0, m);

            // 正交化 (Classical Gram-Schmidt)
            for (int j = 0; j < i; j++) {
                // 计算投影系数
                double coef = doProduct(B[i], B_tilde[j]) / normsSq[j];

                // 对于 m=5304，普通循环利用 CPU 的 SIMD 指令集会比并行流快得多，且不会崩溃。
                for (int k = 0; k < m; k++) {
                    B_tilde[i][k] -= coef * B_tilde[j][k];
                }
            }

            // 更新范数
            normsSq[i] = doProduct(B_tilde[i], B_tilde[i]);
        }

        long endTime = System.currentTimeMillis();
        System.out.println("GSO Time: " + (endTime - startTime) / 1000 + " seconds");
        System.out.println("GSO Completed");
    }

    private static double doProduct(double[] a, double[] b) {
        // Can be parallelized for very large vectors, but overhead usually outweighs benefit for dot product
        double sum = 0;
        for (int i = 0; i < a.length; i++) {
            sum += a[i] * b[i];
        }
        return sum;
    }

    /**
     * Fast SampleD with Precomputed GSO and Parallel Vector Updates.
     * Added `B_double` parameter to avoid repeated doubleValue() conversion inside the loop.
     */
    public static BigInteger[] sampleD_Fast(double[][] B_tilde, double[] normsSq, double[][] B_double, BigIntMatrix B_original, double sigma, double[] c, long nParam) {
        int m = B_original.getRows();
        int n = B_original.getCols();

        BigInteger[] v = new BigInteger[m];
        Arrays.fill(v, BigInteger.ZERO);

        double[] c_current = Arrays.copyOf(c, c.length);

        // Nearest Plane Algorithm
        for (int i = n - 1; i >= 0; i--) {
            double c_prime = doProduct(c_current, B_tilde[i]) / normsSq[i];
            double s_prime = sigma / Math.sqrt(normsSq[i]);

            BigInteger z_i = sampleZ(s_prime, c_prime, nParam);
            double z_double = z_i.doubleValue();

            // Optimization: Parallelize the updates of c_current and v
            final int idx = i;

            // Update c_current: c -= z * b_i
            // Using pre-converted B_double for speed
            IntStream.range(0, m).parallel().forEach(k -> {
                c_current[k] = c_current[k] - z_double * B_double[idx][k];
            });

            // Update v: v += z * b_i
            // This is BigInteger arithmetic, definitely parallelize
            IntStream.range(0, m).parallel().forEach(k -> {
                // Note: BigInteger is immutable, safe for concurrent read of B_original
                // But writing to v array slots is independent, so thread-safe.
                BigInteger val = B_original.get(k, idx).multiply(z_i);
                // v[k] is accessed only by this thread for index k
                // However, v[k] reads previous value.
                // Since we parallelize over 'k', each slot k is independent. Safe.
                v[k] = v[k].add(val);
            });
        }
        return v;
    }

    // Compatibility overload if B_double is not cached
    public static BigInteger[] sampleD_Fast(double[][] B_tilde, double[] normsSq, BigIntMatrix B_original, double sigma, double[] c, long nParam) {
        int n = B_original.getCols();
        int m = B_original.getRows();
        double[][] B_double = new double[n][m];
        IntStream.range(0, n).parallel().forEach(i -> {
            for (int j = 0; j < m; j++) B_double[i][j] = B_original.get(j, i).doubleValue();
        });
        return sampleD_Fast(B_tilde, normsSq, B_double, B_original, sigma, c, nParam);
    }
}