package qut.edu.au.pqcsebaseopt.algorithms;

import qut.edu.au.pqcsebaseopt.tools.BigIntMatrix;
import qut.edu.au.pqcsebaseopt.tools.IncrementalRankChecker;

import java.math.BigInteger;

import static qut.edu.au.pqcsebaseopt.algorithms.SampleZ.sampleZ;

public class SampleR {

    /**
     * SampleR algorithm generate sigma
     * @param m
     * @param q
     * @param sigmaR = sqrt(n * log q)*w(sqrt(log m))
     * @return
     */
    public static BigIntMatrix sampleR(int m, BigInteger q, double sigmaR) {

        BigIntMatrix R = new BigIntMatrix(m, m, q);

        // 优化核心：增量检查器
        IncrementalRankChecker rankChecker = new IncrementalRankChecker(m, q);

        // 逐列生成 (i = 0 到 m-1)
        for (int i = 0; i < m; i++) {

            // 内层循环：直到找到合法的第 i 列
            while (true) {
                BigInteger[] colVector = new BigInteger[m];

                // 采样一列
                for (int j = 0; j < m; j++) {
                    // 调用你的 SampleZ，注意参数顺序
                    colVector[j] = sampleZ(sigmaR, 0.0, m);
                }

                // 增量检查
                if (rankChecker.checkAndAdd(colVector)) {
                    // 合法，填入矩阵
                    for (int j = 0; j < m; j++) {
                        R.set(j, i, colVector[j]);
                    }
                    break; // 进入下一列
                }
                // 否则重试该列
            }
        }

        return R;
    }
}
