package qut.edu.au.pqcsebase.projectsTest.IBEKSTEST;

import ch.qos.logback.core.util.StringCollectionUtil;
import ch.qos.logback.core.util.StringUtil;
import org.junit.Test;
import qut.edu.au.pqcsebase.projects.IBEKS.Setup;
import qut.edu.au.pqcsebase.projects.IBEKS.vo.CT;
import qut.edu.au.pqcsebase.projects.IBEKS.vo.KeyPair;
import qut.edu.au.pqcsebase.projects.IBEKS.vo.MSK;
import qut.edu.au.pqcsebase.projects.IBEKS.vo.TW;
import qut.edu.au.pqcsebase.tools.BigIntMatrix;

import java.sql.SQLOutput;
import java.util.Arrays;

import static qut.edu.au.pqcsebase.projects.IBEKS.Encrypt.encrypt;
import static qut.edu.au.pqcsebase.projects.IBEKS.KeyGen.keyGen;
import static qut.edu.au.pqcsebase.projects.IBEKS.Params.keywords;
import static qut.edu.au.pqcsebase.projects.IBEKS.Test.test;
import static qut.edu.au.pqcsebase.projects.IBEKS.Trapdoor.trapdoor;

public class SetupTest {

    @Test
    public void testSetup() {
        Setup setup = new Setup();
        System.out.println("Setup instance created: " + setup);
    }

    @Test
    public void keyGenTest(){
        Setup setup = new Setup();
        MSK msk = setup.msk;
        String identity = "123456789";
        KeyPair keyPair = keyGen(setup,msk,identity);
        System.out.println("Generated KeyPair for identity " + identity + ": " + keyPair);
        BigIntMatrix result = keyPair.pk_id.multiplyQ(keyPair.sk_id);
        System.out.println("Verification result (A_id * sk_id): " + result);
        System.out.println(result.isZeroMatrix());
    }

    @Test
    public void encryptTest(){
        String keyword = "engineer";
        Setup setup = new Setup();
        MSK msk = setup.msk;
        String identity = "123456789";
        KeyPair keyPair = keyGen(setup,msk,identity);
        //System.out.println("Generated KeyPair for identity " + identity + ": " + keyPair);
        BigIntMatrix result = keyPair.pk_id.multiplyQ(keyPair.sk_id);
        System.out.println("Verification result (A_id * sk_id): " + result);
        System.out.println(result.isZeroMatrix());

        CT ct = encrypt(setup,identity,keyword);
        System.out.println("Generated Ciphertext: " + ct);

    }

    @Test
    public void trapdoorTest(){
        String keyword = "engineer";
        Setup setup = new Setup();
        MSK msk = setup.msk;
        String identity = "123456789";
        KeyPair keyPair = keyGen(setup,msk,identity);
        //System.out.println("Generated KeyPair for identity " + identity + ": " + keyPair);
        BigIntMatrix result = keyPair.pk_id.multiplyQ(keyPair.sk_id);
        System.out.println("Verification result (A_id * sk_id): " + result);
        System.out.println(result.isZeroMatrix());

        CT ct = encrypt(setup,identity,keyword);
        System.out.println("Generated Ciphertext: " + ct);

        //---------------------------Trapdoor-----------------------------------
        String[] keywords = new String[]{"doctor", "engineer", "teacher"};
        TW tw = trapdoor(setup,keywords,keyPair.sk_id,identity);
        System.out.println("Generated Trapdoor: " + tw);

    }

    @Test
    public void testAll(){

        long startTime = System.currentTimeMillis();
        long endTime;

        String keyword = "engineer";
        Setup setup = new Setup();
        MSK msk = setup.msk;
        String identity = "123456789";
        endTime = System.currentTimeMillis();
        System.out.println("-------------------Setup-----------------------------------");
        long setupTime = endTime - startTime;
        System.out.println("Setup Time: " + (endTime - startTime) + " ms");

        startTime = System.currentTimeMillis();
        KeyPair keyPair = keyGen(setup,msk,identity);
        endTime = System.currentTimeMillis();
        System.out.println("-------------------KeyGen-----------------------------------");
        long keyGenTime = endTime - startTime;
        System.out.println("KeyGen Time: " + (endTime - startTime) + " ms");

        //System.out.println("Generated KeyPair for identity " + identity + ": " + keyPair);
        BigIntMatrix result = keyPair.pk_id.multiplyQ(keyPair.sk_id);
        System.out.println("Verification result (A_id * sk_id): " + result);
        System.out.println(result.isZeroMatrix());

        startTime = System.currentTimeMillis();
        CT ct = encrypt(setup,identity,keyword);
        System.out.println("Generated Ciphertext: " + ct);
        endTime = System.currentTimeMillis();
        System.out.println("-------------------Encrypt-----------------------------------");
        long encryptTime = endTime - startTime;
        System.out.println("Encrypt Time: " + (endTime - startTime) + " ms");

        //---------------------------Trapdoor-----------------------------------
        String[] keywords = new String[]{"doctor", "engineer", "teacher"};

        startTime = System.currentTimeMillis();
        TW tw = trapdoor(setup,keywords,keyPair.sk_id,identity);
        endTime = System.currentTimeMillis();
        System.out.println("-------------------Trapdoor-----------------------------------");
        long trapdoorTime = endTime - startTime;
        System.out.println("Trapdoor Time: " + (endTime - startTime) + " ms");
        System.out.println("Generated Trapdoor: " + tw);

        startTime = System.currentTimeMillis();
        Boolean flag = test(ct,tw);
        endTime = System.currentTimeMillis();
        System.out.println("-------------------Test-----------------------------------");
        long testTime = endTime - startTime;
        System.out.println("Test Time: " + (endTime - startTime) + " ms");

        System.out.println("Test result: " + flag);

        System.out.println("----------------------------------------------------------");
        System.out.println("Summary of Times:");
        System.out.println("Setup Time: " + setupTime + " ms");
        System.out.println("KeyGen Time: " + keyGenTime + " ms");
        System.out.println("Encrypt Time: " + encryptTime + " ms");
        System.out.println("Trapdoor Time: " + trapdoorTime + " ms");
        System.out.println("Test Time: " + testTime + " ms");
        System.out.println("----------------------------------------------------------");

    }
}
