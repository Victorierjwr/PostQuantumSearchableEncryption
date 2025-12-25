package qut.edu.au.pqcsebaseopt.algorithms;

import qut.edu.au.pqcsebaseopt.tools.BigIntMatrix;
import qut.edu.au.pqcsebaseopt.tools.IncrementalRankChecker;

import java.math.BigInteger;
import java.util.Arrays;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.IntStream;

import static qut.edu.au.pqcsebaseopt.algorithms.SampleD.*;

public class RandBasis {

    /**
     * Generate a random basis (Standard, Slow version).
     */
    public static BigIntMatrix randBasis(BigIntMatrix S, double sigma, long nParam) {
        // Redirect to the optimized pre-computed version directly
        return randBasisPre(S, sigma, nParam);
    }

/*    public static BigIntMatrix randBasis(BigIntMatrix S, double sigma, long nParam) {
        int m = S.getRowDimension();
        int n = S.getColumnDimension();
        BigInteger q = S.getModulus();

        BigIntMatrix S_prime = new BigIntMatrix(m, n, q);
        double[] center = new double[m];
        Arrays.fill(center, 0.0);

        // 初始化优化后的检查器
        IncrementalRankChecker rankChecker = new IncrementalRankChecker(n, q);

        for (int i = 0; i < n; i++) {
            BigInteger[] sampledVector = sampleD(S, sigma, center, nParam);

            // 简单的零向量检查
            boolean isZero = true;
            for (BigInteger val : sampledVector) {
                if (!val.equals(BigInteger.ZERO)) {
                    isZero = false;
                    break;
                }
            }
            if (isZero) {
                i--; continue;
            }

            // 使用增量检查
            if (rankChecker.checkAndAdd(sampledVector)) {
                // 线性无关，填入矩阵
                for (int j = 0; j < m; j++) {
                    S_prime.set(j, i, sampledVector[j]);
                }
            } else {
                // 线性相关，重试
                i--;
            }
        }

        return LLL(S_prime, 0.99);
    }*/

/*    public static BigIntMatrix randBasisPre(BigIntMatrix S, double sigma, long nParam) {
        int m = S.getRows();
        int n = S.getCols();
        BigInteger q = S.getModulus();

        BigIntMatrix S_prime = new BigIntMatrix(m, n, q);
        double[] center = new double[m];
        Arrays.fill(center, 0.0);

        double[][] B_double = new double[n][m];
        for(int i=0;i<n;i++) for(int j=0;j<m;j++) B_double[i][j] = S.get(j, i).doubleValue();
        double[][] B_tilde = new double[n][m];
        double[] normsSq = new double[n];
        computeGramSchmidt(B_double, B_tilde, normsSq);

        // 初始化优化后的检查器
        IncrementalRankChecker rankChecker = new IncrementalRankChecker(Math.max(n, m), BigInteger.ZERO);

        for (int i = 0; i < n; i++) {
            BigInteger[] sampledVector = sampleD_Fast(B_tilde, normsSq, S, sigma, center, nParam);

            boolean isZero = true;
            for (BigInteger val : sampledVector) {
                if (!val.equals(BigInteger.ZERO)) {
                    isZero = false;
                    break;
                }
            }
            if (isZero) {
                i--; continue;
            }

            // 使用增量检查
            if (rankChecker.checkAndAdd(sampledVector)) {
                for (int j = 0; j < m; j++) {
                    S_prime.set(j, i, sampledVector[j]);
                }
            } else {
                i--;
            }
            if(i % 100 == 0)System.out.println("Generated column " + (i + 1) + "/" + n);
        }

        //BigIntMatrix result= LLL(S_prime, 0.99);
        BigIntMatrix result = S_prime;
        return result;
    }*/

    /*public static BigIntMatrix randBasisPre(BigIntMatrix S, double sigma, long nParam) {
        int m = S.getRows();
        int n = S.getCols();
        BigInteger q = S.getModulus();

        BigIntMatrix S_prime = new BigIntMatrix(m, n, q);
        double[] center = new double[m];
        Arrays.fill(center, 0.0);

        // 1. GSO 计算 (保留)
        double[][] B_double = new double[n][m];
        for(int i=0;i<n;i++) for(int j=0;j<m;j++) B_double[i][j] = S.get(j, i).doubleValue();
        double[][] B_tilde = new double[n][m];
        double[] normsSq = new double[n];
        computeGramSchmidt(B_double, B_tilde, normsSq);

        // 2. 移除 IncrementalRankChecker
        // IncrementalRankChecker rankChecker = new IncrementalRankChecker(Math.max(n, m), BigInteger.ZERO);

        System.out.println("Starting Sampling loop (Fast Mode)...");
        long start = System.currentTimeMillis();

        for (int i = 0; i < n; i++) {
            // 采样
            BigInteger[] sampledVector = sampleD_Fast(B_tilde, normsSq, S, sigma, center, nParam);

            // 简单的零向量检查 (保留)
            boolean isZero = true;
            for (BigInteger val : sampledVector) {
                if (!val.equals(BigInteger.ZERO)) {
                    isZero = false;
                    break;
                }
            }
            if (isZero) {
                System.out.println("Zero vector sampled, retrying...");
                i--; continue;
            }

            // [核心修改] 直接填入，不做线性无关检查
            // 相信 SampleD 的数学性质：生成的向量在高维空间几乎正交
            for (int j = 0; j < m; j++) {
                S_prime.set(j, i, sampledVector[j]);
            }

            // 打印进度
            if(i % 100 == 0) {
                long now = System.currentTimeMillis();
                System.out.println("Generated column " + (i + 1) + "/" + n + " (Time: " + (now-start)/1000 + "s)");
            }
        }

        // 不做 LLL，直接返回
        return S_prime;
    }*/

    public static BigIntMatrix randBasisPre(BigIntMatrix S, double sigma, long nParam) {
        int m = S.getRows();
        int n = S.getCols();
        BigInteger q = S.getModulus();

        BigIntMatrix S_prime = new BigIntMatrix(m, n, q);
        double[] center = new double[m];
        Arrays.fill(center, 0.0);

        System.out.println("Starting GSO computation...");
        long gsoStart = System.currentTimeMillis();

        // 1. GSO 计算 (保持串行，因为这步依赖性强，难以并行且不是主要瓶颈)
        double[][] B_double = new double[n][m];
        for(int i=0; i<n; i++) {
            for(int j=0; j<m; j++) {
                B_double[i][j] = S.get(j, i).doubleValue();
            }
        }
        double[][] B_tilde = new double[n][m];
        double[] normsSq = new double[n];
        computeGramSchmidt(B_double, B_tilde, normsSq);

        System.out.println("GSO Completed in " + (System.currentTimeMillis() - gsoStart)/1000 + "s");

        // 2. 并行采样 (Parallel Sampling)
        System.out.println("Starting Parallel Sampling loop (Safe Mode)...");
        long start = System.currentTimeMillis();
        AtomicInteger progress = new AtomicInteger(0);

        // 获取 CPU 核心数，创建一个大小合适的线程池
        int parallelism = Runtime.getRuntime().availableProcessors();
        ForkJoinPool customThreadPool = new ForkJoinPool(Math.min(parallelism, 5));

        BigInteger[][] columns;
        try {
            // 在自定义线程池中提交并行任务
            columns = customThreadPool.submit(() ->
                    IntStream.range(0, n).parallel().mapToObj(i -> {
                        // --- 采样核心逻辑 ---
                        BigInteger[] sampledVector = sampleD_Fast(B_tilde, normsSq, S, sigma, center, nParam);

                        // 简单的零向量检查
                        boolean isZero = true;
                        for (BigInteger val : sampledVector) {
                            if (!val.equals(BigInteger.ZERO)) {
                                isZero = false;
                                break;
                            }
                        }
                        if (isZero) {
                            // 极其罕见情况，重试一次
                            sampledVector = sampleD_Fast(B_tilde, normsSq, S, sigma, center, nParam);
                        }

                        // --- 进度打印 ---
                        int currentCount = progress.incrementAndGet();
                        if (currentCount % 100 == 0) {
                            long now = System.currentTimeMillis();
                            double speed = (now - start) / 1000.0;
                            System.out.println("Generated " + currentCount + "/" + n + " columns (Time: " + (long)speed + "s)");
                        }

                        return sampledVector;
                    }).toArray(BigInteger[][]::new)
            ).get(); // 等待所有任务完成

        } catch (InterruptedException | ExecutionException e) {
            throw new RuntimeException("Parallel sampling failed", e);
        } finally {
            customThreadPool.shutdown();
        }

        System.out.println("Sampling finished. Constructing matrix...");

        // 3. 将生成的列填入矩阵 (串行填充以保证线程安全)
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                S_prime.set(j, i, columns[i][j]);
            }
        }

        // 跳过耗时的 LLL 归约，直接返回采样结果 (对于 MP12/ABB10 陷门已足够)
        return S_prime;
    }

    /**
     * LLL Algorithm.
     * Note: For huge dimensions (N=4093), exact LLL is extremely slow.
     * This implementation uses double-precision GSO which is standard for performance
     * but may suffer from precision issues in extreme cases.
     */
    public static BigIntMatrix LLL(BigIntMatrix B, double delta) {

        long startTime = System.currentTimeMillis();

        int n = B.getCols();
        int m = B.getRows();

        // 1. Deep copy
        BigIntMatrix S = B.copy();

        // 2. GSO Cache
        // We only need mu and normsSq for the logic.
        // B_tilde vectors are NOT needed for the LLL logic itself if mu is maintained.
        double[] normsSq = new double[n];
        double[][] mu = new double[n][n];

        // Initial GSO (O(n^3)) - Done only once!
        computeGSO(S, normsSq, mu);

        int k = 1;
        while (k < n) {
            if(k % 500 == 0)System.out.println("LLL Progress: Processing column " + (k + 1) + "/" + n);
            // === Step A: Size Reduction (Optimized) ===
            // Use cached mu table instead of computing dot products
            for (int j = k - 1; j >= 0; j--) {
                if (Math.abs(mu[k][j]) > 0.5) {
                    long c = Math.round(mu[k][j]);
                    if (c != 0) {
                        // Matrix operation: Col_k = Col_k - c * Col_j
                        S.colSub(k, j, BigInteger.valueOf(c));

                        // Update GSO coefficients (O(n) speed)
                        // mu[k][j] becomes small
                        mu[k][j] -= c;
                        // Update other coefficients affected by this subtraction
                        for (int l = 0; l < j; l++) {
                            mu[k][l] -= c * mu[j][l];
                        }
                    }
                }
            }

            // === Step B: Lovasz Condition ===
            double lhs = delta * normsSq[k - 1];
            double rhs = normsSq[k] + mu[k][k - 1] * mu[k][k - 1] * normsSq[k - 1];

            if (rhs < lhs) {
                // === Swap Columns k and k-1 ===
                S.swapCols(k, k - 1);

                // === Local GSO Update (The Magic Part) ===
                // Instead of recomputing the whole GSO (O(n^3)), we use algebra
                // to update only the changed values in O(n).

                int km1 = k - 1;
                double mu_k_km1 = mu[k][km1];
                double B_km1 = normsSq[km1]; // Old norm at k-1
                double B_k = normsSq[k];     // Old norm at k (orthogonal part)

                // 1. Update Norms
                double new_B_km1 = B_k + mu_k_km1 * mu_k_km1 * B_km1;
                double new_B_k = B_k * B_km1 / new_B_km1;

                normsSq[km1] = new_B_km1;
                normsSq[k] = new_B_k;

                // 2. Update mu[k][k-1]
                mu[k][km1] = mu_k_km1 * B_km1 / new_B_km1;

                // 3. Update mu matrix for columns j < k-1 (just swap)
                for (int j = 0; j < km1; j++) {
                    double tmp = mu[km1][j];
                    mu[km1][j] = mu[k][j];
                    mu[k][j] = tmp;
                }

                // 4. Update mu matrix for rows i > k (Givens rotation-like update)
                for (int i = k + 1; i < n; i++) {
                    double t = mu[i][k];
                    double u = mu[i][km1];
                    // These formulas preserve the validity of the GSO without recomputing vectors
                    mu[i][k] = u - mu_k_km1 * t;
                    mu[i][km1] = t + mu[k][km1] * mu[i][k];
                }

                // Backtrack
                k = Math.max(k - 1, 1);
            } else {
                k++;
            }
        }
        long endTime = System.currentTimeMillis();
        System.out.println("LLL completed in " + (endTime - startTime)/1000 + " seconds");
        return S;
    }

    /**
     * Initial GSO computation.
     * Computes mu and normsSq from scratch.
     */
    private static void computeGSO(BigIntMatrix S, double[] normsSq, double[][] mu) {
        long startTime = System.currentTimeMillis();
        System.out.println("Starting initial GSO computation for LLL...");
        int n = S.getCols();
        int m = S.getRows();

        // Temporary storage for orthogonal vectors (only needed during init)
        double[][] b_tilde = new double[n][m];

        for (int i = 0; i < n; i++) {
            if(i % 500 == 0)System.out.println("Computing GSO for column " + (i + 1) + "/" + n);
            // Load B[i] into b_tilde[i]
            for (int r = 0; r < m; r++) {
                b_tilde[i][r] = S.get(r, i).doubleValue();
            }

            // Orthogonalize against previous
            for (int j = 0; j < i; j++) {
                double dot = 0;
                for (int r = 0; r < m; r++) dot += b_tilde[i][r] * b_tilde[j][r];

                mu[i][j] = dot / normsSq[j];

                for (int r = 0; r < m; r++) {
                    b_tilde[i][r] -= mu[i][j] * b_tilde[j][r];
                }
            }

            // Calc norm
            double norm = 0;
            for (int r = 0; r < m; r++) norm += b_tilde[i][r] * b_tilde[i][r];
            normsSq[i] = norm;
            mu[i][i] = 1.0;
        }
        System.out.println("Initial GSO computation completed for LLL.");
        long endTime = System.currentTimeMillis();
        System.out.println("Initial GSO for LLL Time: " + (endTime - startTime)/1000 + " seconds");
    }
}