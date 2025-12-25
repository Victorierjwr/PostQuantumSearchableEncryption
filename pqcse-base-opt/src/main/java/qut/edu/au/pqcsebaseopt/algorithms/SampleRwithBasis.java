package qut.edu.au.pqcsebaseopt.algorithms;


import qut.edu.au.pqcsebaseopt.tools.BigIntMatrix;
import qut.edu.au.pqcsebaseopt.tools.IncrementalRankChecker;

import java.math.BigInteger;

import static qut.edu.au.pqcsebaseopt.algorithms.SamplePre.samplePre;
import static qut.edu.au.pqcsebaseopt.algorithms.TrapGen.trapGen;


public class SampleRwithBasis {

    /**
     * Algorithm SampleRwithBasis(A) from paper.
     * 1. Generate (B, Tb) using TrapGen.
     * 2. Sample columns r_i using SamplePre such that B*R = A mod q.
     * 3. Ensure columns of R are linearly independent using IncrementalRankChecker.
     */
    public static BigIntMatrix[] sampleRwithBasis(BigIntMatrix A, double sigmaR) {
        int n = A.getRows();
        int m = A.getCols();
        BigInteger q = A.getModulus();

        // === Step 1: Run TrapGen(q, n) ===
        BigIntMatrix[] trapdoor = trapGen(n, q.longValue(), 1);
        BigIntMatrix B = trapdoor[0];
        BigIntMatrix Tb = trapdoor[1];

        BigIntMatrix R = new BigIntMatrix(m, m, q);

        //check rank
        IncrementalRankChecker rankChecker = new IncrementalRankChecker(m, q);

        // === Step 2: Loop i = 1 to m ===
        for (int i = 0; i < m; i++) {
            // 获取目标向量 u = A的第i列
            BigInteger[] u = new BigInteger[n];
            for (int k = 0; k < n; k++) u[k] = A.get(k, i);

            // (2b) Repeat until linearly independent
            while (true) {
                // (2a) sample r_i using SamplePre
                BigInteger[] r_i = samplePre(B, Tb, sigmaR, u);

                // 检查线性无关性
                if (rankChecker.checkAndAdd(r_i)) {
                    // 成功：填入矩阵 R
                    for (int k = 0; k < m; k++) {
                        R.set(k, i, r_i[k]);
                    }
                    break; // 跳出 while，处理下一列
                }
                // 失败：线性相关，自动重试 (while 循环继续)
                // System.out.println("Dependent vector at col " + i + ", retrying...");
            }
        }

        // === Step 3: Output R and Tb ===
        // 为了方便验证 B*R=A，通常也会返回 B
        return new BigIntMatrix[]{R, Tb, B};
    }
}