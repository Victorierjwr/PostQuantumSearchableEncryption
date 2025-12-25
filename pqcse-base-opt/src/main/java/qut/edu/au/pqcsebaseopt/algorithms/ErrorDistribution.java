package qut.edu.au.pqcsebaseopt.algorithms;

import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.MathContext;
import java.util.concurrent.ThreadLocalRandom;

public class ErrorDistribution {

    // 使用密码学安全的随机数生成器
    private static final ThreadLocalRandom secureRandom = ThreadLocalRandom.current();

    private static final double PI = 3.14159;

    /**
     * 生成服从高斯分布的离散误差项，并映射到 Z_q
     *
     * @param q 模数 (对应 NTL::ZZ)
     */
    public static BigInteger generateError(BigInteger q) {
        double a = 0.05; // 基础参数 (类似 LWE 的 alpha)

        // 计算方差/缩放因子 D
        // C++: D = (a / sqrt(2 * PI))*(a / sqrt(2 * PI));
        double term = a / Math.sqrt(2 * PI);

        double D = term * term;  //source code
        //double D = term;  //recommend

        // 2. 生成均匀分布随机数 u1, u2
        double u1 = secureRandom.nextDouble();
        double u2 = secureRandom.nextDouble();

        // 避免 log(0) 的情况
        while (u1 <= 0.0) {
            u1 = secureRandom.nextDouble();
        }

        // 3. Box-Muller 变换生成标准正态分布变量 C
        // A = sqrt((-2)*log(uni[0]));
        double A = Math.sqrt(-2.0 * Math.log(u1));
        // B = 2 * PI*uni[1];
        double B = 2.0 * PI * u2;
        // C = A * cos(B); (标准正态分布 N(0,1))
        double C = A * Math.cos(B);

        // 4. 缩放
        // r1 = E + C * D; (C++中 E=0)
        double r1 = C * D;

        // 检查无穷大
        if (Double.isInfinite(r1) || Double.isNaN(r1)) {
            // 极低概率事件，递归重试或返回 0
            return generateError(q);
        }

        // 5. 映射到整数并取模
        // 逻辑: r2 = round(q * r1) mod q

        // 将 double r1 转为 BigDecimal 以保持与 BigInteger q 运算的精度
        BigDecimal r1BD = new BigDecimal(r1);
        BigDecimal qBD = new BigDecimal(q);

        // 计算 q * r1
        BigDecimal product = qBD.multiply(r1BD);

        // 四舍五入 (round)
        BigInteger roundedValue = product.round(MathContext.DECIMAL128).toBigInteger(); // 默认是四舍五入

        // 取模 (mod q)
        return roundedValue.mod(q);
    }

    /**
     * 测试主函数
     */
    public static void main(String[] args) {
        BigInteger q = new BigInteger("10007"); // 示例模数
        System.out.println("Generating errors mod " + q + ":");

        for (int i = 0; i < 10; i++) {
            System.out.println(generateError(q));
        }
    }
}