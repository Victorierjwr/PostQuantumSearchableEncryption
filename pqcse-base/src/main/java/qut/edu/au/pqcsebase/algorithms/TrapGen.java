package qut.edu.au.pqcsebase.algorithms;

import qut.edu.au.pqcsebase.tools.BigIntMatrix;

import java.math.BigInteger;

public class TrapGen {

    /**
     * flag = 1 : first method
     * flag = 2 : second method
     * Main trapdoor generation function.
     * Returns [A1, S] where A1 is n x m1 matrix and S is (m1 + m2) x (m1 + m2) matrix.
     */
    public static BigIntMatrix[] trapGen(long n, long q, int flag){
        BigIntMatrix[] A_S = new BigIntMatrix[2];
        if(flag == 1){
            A_S = trapGenFirst(n,q);
        }
        if(flag == 2){
            A_S = trapGenSecond(n,q);
        }
        return A_S;
    }

    private static BigIntMatrix[] trapGenFirst(long n, long q) {

        BigIntMatrix[] A_S = new BigIntMatrix[2];

        TrapGenFirst TGF = new TrapGenFirst(n, q);
        int m1 = TGF.getM1();
        int m2 = TGF.getM2();
        int m = TGF.getM();
        int l = TGF.getL();
        BigInteger qBig = BigInteger.valueOf(q);
        // 1. Generate A1
        BigIntMatrix A1 = TGF.genA1();

        // 2. Compute Basis / HNF
        // For testing purposes, we construct a full rank matrix to get H
        BigIntMatrix X = BigIntMatrix.image(A1);
        BigIntMatrix M = X.getSubMatrix(0, (int)n-1, 0, (int)n-1);

        //M: n*n A1_M: n*(m1 - n) _A1_M: n*(m1 - n)
        BigIntMatrix O = BigIntMatrix.zeroMatrix((int)n, m1 - (int)n, qBig);
        BigIntMatrix A1_M = X.getSubMatrix(0, (int)n-1, (int) n, m1 - 1);

        // -A1_M
        BigIntMatrix _A1_M = O.subtract(A1_M);
//        TGF.print("-A1_M", _A1_M);

        BigIntMatrix M_1 = BigIntMatrix.inverseQ(M);
//        TGF.print("M^-1", M_1);

        // x1 = n*(m1 - n)
        BigIntMatrix x1 = M_1.multiplyQ(_A1_M);
//        TGF.print("x1", x1);

        // x2 = (m1 - n)*(m1 - n)
        BigIntMatrix x2 = BigIntMatrix.identity(m1 - (int)n, qBig);
//        TGF.print("x2", x2);

        // x1_x2 = m1 * (m1 - n)
        BigIntMatrix x1_x2 = BigIntMatrix.rowConcat(x1, x2);
//        TGF.print("x1_x2", x1_x2);

        // x3 = m1 *  m1  need get m1 * n
        BigIntMatrix x3 = BigIntMatrix.diagonal(m1, qBig, qBig);
//        TGF.print("x3", x3);

        BigIntMatrix beforeH = BigIntMatrix.columnConcat(x1_x2,x3.getSubMatrix(0, m1 - 1, m1 - (int)n , m1 - 1));
//       TGF.print("beforeH", beforeH);

        BigIntMatrix H = TGF.genH(beforeH); // Returns BigIntMatrix, no overflow!
//        TGF.print("H", H);


        //--------------------------Construction------------------------------------
        //--------------------------Construction C----------------------------------
        BigIntMatrix I = BigIntMatrix.identity(m1, qBig);

        BigIntMatrix C = I;
//        TGF.print("C", C);

        //--------------------------Construction G----------------------------------
        // test -------  AH = A1 * H
        BigIntMatrix AH = A1.multiply(H);
//        TGF.print("AH", AH);

        // H‘ = H - I
        BigIntMatrix H_I = H.subtractQ(I); // H - I
//        TGF.print("H - I", H_I);

        // G m1*m2
        BigIntMatrix G = TGF.genG(H_I);
//        TGF.print("G", G);

        //--------------------------Construction P----------------------------------
        BigIntMatrix P = TGF.genP();
//        TGF.print("P", P);

        //verify GP = H_I
        BigIntMatrix GP = G.multiply(P);
//        TGF.print("GP", GP);
//        System.out.println("GP == H_I: " + GP.equals(H_I));

        //--------------------------Construction U----------------------------------
        BigIntMatrix Tl = TGF.genTl(l, qBig);
//        TGF.print("Tl", Tl);

        BigIntMatrix U = TGF.genU();
//        TGF.print("U", U);
        //--------------------------Construction R----------------------------------
        BigIntMatrix R = TGF.genR();
//        TGF.print("R", R);

        //Verify A1(GP + C) = 0 mod q  n*m1
        BigIntMatrix GP_C = GP.add(C);
        BigIntMatrix AGP_C = A1.multiplyQ(GP_C);
//        TGF.print("A1(GP + C)", AGP_C);
//        System.out.println("A1(GP + C) == 0 mod q: " + AGP_C.isZeroMatrix());

        //--------------------------Algorithm Part Output A S----------------------------------

        //--------------------------Construction A = [A1 | A2]----------------------------------
        BigIntMatrix RG = R.add(G);
        BigIntMatrix tmpA2 = A1.multiplyQ(RG);
        BigIntMatrix tmp0 = BigIntMatrix.zeroMatrix((int)n, m2, qBig);
        BigIntMatrix A2 = tmp0.subtractQ(tmpA2);
//        TGF.print("A2", A2);

        BigIntMatrix A = BigIntMatrix.columnConcat(A1, A2);
        A_S[0] = A;
//        TGF.print("A", A);


        //--------------------------Construction S four part (G+R)U RP - C U P  no mod q----------------------------------
        //1.left_up  (G+R)U
        BigIntMatrix GR_U = (G.add(R)).multiply(U);
//        TGF.print("(G+R)U", GR_U);
        //2.right_up  RP - C
        BigIntMatrix RP_C = (R.multiply(P)).subtract(C);
//        TGF.print("RP - C", RP_C);
        //3.left_down U
        BigIntMatrix U_down = U;
//        TGF.print("U", U_down);
        //4.right_down P
        BigIntMatrix P_down = P;
//        TGF.print("P", P_down);

        //---------Combination S = [ (G+R)U | RP - C ; U | P ] (m1 + m2)) * (m1 + m2)
        BigIntMatrix upper = BigIntMatrix.columnConcat(GR_U, RP_C);
        BigIntMatrix lower = BigIntMatrix.columnConcat(U_down, P_down);
        BigIntMatrix S = BigIntMatrix.rowConcat(upper, lower);
        A_S[1] = S;
//        TGF.print("S", S);
        return A_S;
    }

    private static BigIntMatrix[] trapGenSecond(long n, long q) {
        return null;
    }
}
