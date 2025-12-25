package qut.edu.au.pqcsebase.algorithmTest;

import org.junit.Test;
import qut.edu.au.pqcsebase.algorithms.BasisDel;
import qut.edu.au.pqcsebase.algorithms.TrapGen;
import qut.edu.au.pqcsebase.tools.BigIntMatrix;

import java.math.BigInteger;

import static java.lang.Math.log;
import static java.lang.Math.sqrt;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static qut.edu.au.pqcsebase.algorithms.SampleR.sampleR;
import static qut.edu.au.pqcsebase.algorithms.BasisDel.basisDel;
import static qut.edu.au.pqcsebase.algorithms.TrapGen.trapGen;
import static qut.edu.au.pqcsebase.tools.BigIntMatrix.*;


public class BasisDelTest {

    @Test
    public void testSampleRProperties() {
        System.out.println("=== Starting SampleR Test ===");

        // 1. 设置参数
        int m = 50;              // 维度 (建议测试 20~100)
        BigInteger q = BigInteger.valueOf(408); // 大素数模数
        double sigmaR = 4.0;     // 高斯参数 (通常较小)

        System.out.println("Parameters: m=" + m + ", q=" + q + ", sigmaR=" + sigmaR);

        // 2. 运行 SampleR
        long startTime = System.currentTimeMillis();
        BigIntMatrix R = sampleR(m, q, sigmaR);
        long endTime = System.currentTimeMillis();

        System.out.println("SampleR finished in " + (endTime - startTime) + " ms");

        // === 验证 1: 维度检查 ===
        assertEquals("Rows must be m", m, R.getRowDimension());
        assertEquals("Cols must be m", m, R.getColumnDimension());

        // === 验证 2: 模 q 下的可逆性 (满秩) ===
        // 这是 SampleR 最核心的要求
        System.out.println("Checking Rank mod q...");
        int rank = gaussMod(R,q);
        System.out.println("Rank: " + rank + " (Expected: " + m + ")");
        assertEquals("Matrix R must be full rank (invertible) mod q", m, rank);

        // === 验证 3: 低范数 (元素大小检查) ===
        // 检查元素是否都在合理的范围内 (例如 +/- 6*sigma)
        // 对于 sigma=4.0, 绝大多数元素应该在 [-24, 24] 之间
        System.out.println("Checking Entry Norms...");

        BigInteger maxVal = BigInteger.ZERO;
        BigInteger minVal = BigInteger.ZERO;
        double bound = sigmaR * 6.0;
        boolean isLowNorm = true;

        // 统计分布情况
        int zeroCount = 0;

        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                BigInteger val = R.get(i, j);

                if (val.compareTo(maxVal) > 0) maxVal = val;
                if (val.compareTo(minVal) < 0) minVal = val;

                if (val.equals(BigInteger.ZERO)) zeroCount++;

                // 检查是否越界 (虽然高斯分布有极小概率越界，但测试中一般不应出现极端值)
                if (Math.abs(val.doubleValue()) > bound) {
                    System.err.println("Warning: Outlier found at (" + i + "," + j + "): " + val);
                    // 注意：这里用 warn 而不是 fail，因为统计学上极值是可能的
                }
            }
        }

        System.out.println("Value Range: [" + minVal + ", " + maxVal + "]");
        System.out.println("Zero Frequency: " + zeroCount + "/" + (m*m) +
                " (" + String.format("%.2f", (double)zeroCount/(m*m)*100) + "%)");

        // 验证范围是否合理 (对于 sigma=4，最大值通常不会超过 30)
        assertTrue("Values should be small integers", maxVal.intValue() < 50 && minVal.intValue() > -50);

        System.out.println("=== Test Passed Successfully ===");
    }

    @Test
    public void basisDelTest() {

        long n = 5;
        long q = 17;
        int flag = 1;
        BigIntMatrix[] A_S = trapGen(n,q,flag);
        BigIntMatrix A = A_S[0];
        BigIntMatrix S = A_S[1];
        int m = S.getRowDimension();
        System.out.println("m :" + m);
        double sigmaR = sqrt(n * log(q)) * sqrt(log(m));
        double sigma = 10000;
        System.out.println("sigmaR: " + sigmaR);
        BigIntMatrix R = sampleR(m, BigInteger.valueOf(q),sigmaR);
        System.out.println("R: " + R.toString());

        System.out.println("rankMod R: " + gaussMod(R,BigInteger.valueOf(q)));
        System.out.println("rank R: " + gaussMod(R,BigInteger.ZERO));

        BigIntMatrix B = A.multiplyQ(inverseQ(R));
        System.out.println("B: " + B.toString());

        BigIntMatrix TB = basisDel(A, R, S, sigma, n);
        System.out.println("TB: " + TB.toString());

        System.out.println("Verify B * TB = 0 mod q: " );
        BigIntMatrix B_TB = B.multiplyQ(TB);
        System.out.println("B_TB: " + B_TB.toString());
        System.out.println("B * TB == 0 mod q: " + B_TB.isZeroMatrix());

    }

}
