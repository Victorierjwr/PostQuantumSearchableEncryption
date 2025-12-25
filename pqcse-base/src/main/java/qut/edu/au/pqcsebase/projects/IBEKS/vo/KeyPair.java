package qut.edu.au.pqcsebase.projects.IBEKS.vo;

import lombok.*;
import qut.edu.au.pqcsebase.tools.BigIntMatrix;

@Data
@Setter @Getter
@AllArgsConstructor@NoArgsConstructor
@ToString
public class KeyPair {

    public BigIntMatrix pk_id;
    public BigIntMatrix sk_id;

}
