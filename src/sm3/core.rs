// 初始值，用于确定压缩函数寄存器的初态
const IV: [u32; 8] = [0x7380_166f, 0x4914_b2b9, 0x1724_42d7, 0xda8a_0600, 0xa96f_30bc, 0x1631_38aa, 0xe38d_ee4d, 0xb0fb_0e4e];

const T0: u32 = 0x79cc_4519;
const T1: u32 = 0x7a87_9d8a;

fn ff0(x: u32, y: u32, z: u32) -> u32 {
    x ^ y ^ z
}

fn ff1(x: u32, y: u32, z: u32) -> u32 {
    (x & y) | (x & z) | (y & z)
}

fn gg0(x: u32, y: u32, z: u32) -> u32 {
    x ^ y ^ z
}

fn gg1(x: u32, y: u32, z: u32) -> u32 {
    (x & y) | (!x & z)
}

/// 压缩函数中的置换函数
fn p0(x: u32) -> u32 {
    x ^ x.rotate_left(9) ^ x.rotate_left(17)
}

/// 消息扩展中的置换函数
fn p1(x: u32) -> u32 {
    x ^ x.rotate_left(15) ^ x.rotate_left(23)
}


#[derive(Debug)]
pub struct Crypto {
    data: Vec<u8>,
    blocks: Vec<[u8; 64]>,
    registers: [u32; 8],
}

impl Crypto {
    pub fn new(data: &[u8]) -> Self {
        Crypto {
            data: data.iter().map(|e| *e).collect(),
            blocks: Vec::new(),
            registers: IV,
        }
    }

    pub fn hash(&mut self) -> [u8; 32] {
        self.pad().block().iterate().output()
    }

    /// 假设消息m的长度为l 比特。首先将比特“1”添加到消息的末尾，再添加k 个“0”，
    /// k是满足l + 1 + k ≡ 448mod512 的最小的非负整数。然后再添加一个64位比特串，该比特串是长度l的二进 制表示。
    /// 填充后的消息m′的比特长度为512的倍数。
    /// 例如:对消息01100001 01100010 01100011，其长度l=24，经填充得到比特串:
    /// 01100001 01100010 01100011 1 {00 · · · 00}(423比特) {00 · · · 011000}(64比特，l的二进制表示)
    fn pad(&mut self) -> &mut Self {
        // 计算原始数据的长度
        let l = (self.data.len() << 3) as u64;
        // 将'10000000'添加到数据的末尾
        self.data.push(0x80);
        // 循环n次填充0x00, l + 8 + k = 448 mod 512,  k mod 8 = n
        while self.data.len() % 64 != 56 {
            self.data.push(0x00);
        }
        // 填充l的二进制表示，长度64位；填充后的数据总长度为512 * N位。
        self.data.push((l >> 56 & 0xff) as u8);
        self.data.push((l >> 48 & 0xff) as u8);
        self.data.push((l >> 40 & 0xff) as u8);
        self.data.push((l >> 32 & 0xff) as u8);
        self.data.push((l >> 24 & 0xff) as u8);
        self.data.push((l >> 16 & 0xff) as u8);
        self.data.push((l >> 8 & 0xff) as u8);
        self.data.push((l & 0xff) as u8);
        self
    }


    /// 分组： 将填充后的消息m′按512比特进行分组:m′ = B(0)B(1) · · · B(n−1), 其中n=(l+k+65)/512。
    fn block(&mut self) -> &mut Self {
        let length = self.data.len();
        let mut c = 0;
        while c * 64 != length {
            let mut block = [0; 64];
            for i in (c * 64)..((c + 1) * 64) {
                block[i - c * 64] = self.data[i];
            }
            self.blocks.push(block);
            c += 1;
        }
        self
    }

    /// 迭代压缩
    /// 1. 扩展
    ///     将消息分组B(i)按以下方法扩展生成132个字W0, W1, · · · , W67, W0′, W1′, · · · , W63′，
    ///     用于压缩函数CF:
    ///         a)将消息分组B(i)划分为16个字W0, W1, · · · , W15。
    ///         b)FOR j=16 TO 67
    ///             Wj ← P1(Wj−16 ⊕Wj−9 ⊕(Wj−3 ≪ 15))⊕(Wj−13 ≪ 7)⊕Wj−6
    ///         c)FOR j=0 TO 63
    ///              Wj′ =Wj ⊕Wj+4
    /// 2.压缩
    ///     令A,B,C,D,E,F,G,H为字寄存器,SS1,SS2,TT1,TT2为中间变量,压缩函数V i+1 i ≤ n−1。
    ///     计算过程描述如下:
    ///     ABCDEFGH ← V (i)
    ///     FOR j=0 TO 63
    ///         SS1←((A≪12)+E+(Tj ≪j))≪7 SS2 ← SS1⊕(A ≪ 12)
    ///         TT1 ← FFj (A, B, C) + D + SS2 + Wj′
    ///         TT2 ← GGj (E, F, G) + H + SS1 + Wj
    ///         D←C
    ///         C←B≪9
    ///         B←A
    ///         A←TT1
    ///         H←G
    ///         G ← F ≪ 19
    ///         F←E
    ///         E ← P0(TT2)
    ///     V(i+1) ← ABCDEFGH⊕V(i)
    fn iterate(&mut self) -> &mut Self {
        self.blocks.iter().for_each(|b| {
            // 扩展
            // 每个分组扩展生成132个字W0, W1, · · · , W67, W0′, W1′, · · · , W63′
            let mut w1: [u32; 68] = [0; 68];
            let mut w2: [u32; 64] = [0; 64];
            // 将消息分组B(i)划分为16个字 W0, W1, · · · , W15
            for i in 0..16 {
                w1[i] = u32::from(b[i * 4]) << 24
                    | u32::from(b[i * 4 + 1]) << 16
                    | u32::from(b[i * 4 + 2]) << 8
                    | u32::from(b[i * 4 + 3]);
            }
            // 计算 W16, ..., W67;  Wj ← P1(Wj−16 ⊕ Wj−9 ⊕ (Wj−3 ≪ 15)) ⊕ (Wj−13 ≪ 7) ⊕ Wj−6
            for i in 16..68 {
                w1[i] = p1(w1[i - 16] ^ w1[i - 9] ^ w1[i - 3].rotate_left(15))
                    ^ w1[i - 13].rotate_left(7)
                    ^ w1[i - 6];
            }
            // 计算 W': W'0, W'1, ... W'63;   Wj′ = Wj ⊕ Wj+4
            for i in 0..64 {
                w2[i] = w1[i] ^ w1[i + 4];
            }
            // 压缩
            // ABCDEFGH ← V (i)
            let mut ra = self.registers[0];
            let mut rb = self.registers[1];
            let mut rc = self.registers[2];
            let mut rd = self.registers[3];
            let mut re = self.registers[4];
            let mut rf = self.registers[5];
            let mut rg = self.registers[6];
            let mut rh = self.registers[7];

            let mut ss1: u32;
            let mut ss2: u32;
            let mut tt1: u32;
            let mut tt2: u32;
            for i in 0..16 {
                ss1 = ra.rotate_left(12)
                    .wrapping_add(re)
                    .wrapping_add(T0.rotate_left(i as u32))
                    .rotate_left(7);
                ss2 = ss1 ^ ra.rotate_left(12);
                tt1 = ff0(ra, rb, rc)
                    .wrapping_add(rd)
                    .wrapping_add(ss2)
                    .wrapping_add(w2[i]);
                tt2 = gg0(re, rf, rg)
                    .wrapping_add(rh)
                    .wrapping_add(ss1)
                    .wrapping_add(w1[i]);
                rd = rc;
                rc = rb.rotate_left(9);
                rb = ra;
                ra = tt1;
                rh = rg;
                rg = rf.rotate_left(19);
                rf = re;
                re = p0(tt2);
            }
            for i in 16..64 {
                ss1 = ra.rotate_left(12)
                    .wrapping_add(re)
                    .wrapping_add(T1.rotate_left(i as u32))
                    .rotate_left(7);
                ss2 = ss1 ^ ra.rotate_left(12);
                tt1 = ff1(ra, rb, rc)
                    .wrapping_add(rd)
                    .wrapping_add(ss2)
                    .wrapping_add(w2[i]);
                tt2 = gg1(re, rf, rg)
                    .wrapping_add(rh)
                    .wrapping_add(ss1)
                    .wrapping_add(w1[i]);
                rd = rc;
                rc = rb.rotate_left(9);
                rb = ra;
                ra = tt1;
                rh = rg;
                rg = rf.rotate_left(19);
                rf = re;
                re = p0(tt2);
            }
            // V(i+1) ← ABCDEFGH⊕V(i)
            self.registers[0] ^= ra;
            self.registers[1] ^= rb;
            self.registers[2] ^= rc;
            self.registers[3] ^= rd;
            self.registers[4] ^= re;
            self.registers[5] ^= rf;
            self.registers[6] ^= rg;
            self.registers[7] ^= rh;
        });
        self
    }

    /// 输出256比特的哈希值
    fn output(&self) -> [u8; 32] {
        // 大端模式：[u32; 8] -> [u8; 32]
        let mut hash: [u8; 32] = [0; 32];
        for (i, e) in self.registers.iter().enumerate() {
            hash[i * 4] = (*e >> 24) as u8;
            hash[i * 4 + 1] = (*e >> 16) as u8;
            hash[i * 4 + 2] = (*e >> 8) as u8;
            hash[i * 4 + 3] = *e as u8;
        }
        hash
    }
}


#[cfg(test)]
mod tests {
    use crate::sm3::core::Crypto;

    #[test]
    fn main() {
        let plain = String::from("abc");
        let data = plain.as_bytes();
        let hash = hex::encode(Crypto::new(data).hash());
        assert_eq!(hash, "66c7f0f462eeedd9d1f2d46bdc10e4e24167c4875cf2f7a2297da02b8f4ba8e0");
    }
}


