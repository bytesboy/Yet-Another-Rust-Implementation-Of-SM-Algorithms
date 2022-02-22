pub mod sm2;
pub mod sm3;
pub mod sm4;


#[cfg(test)]
mod tests {
    use crate::{sm2, sm3, sm4};


    #[test]
    fn sm2_keypair() {
        for i in 0..10 {
            let (prk, puk) = sm2::generate_key();
            println!("{:?}: ({:?}, {:?})", i, prk, puk);
        }
    }

    #[test]
    fn sm2_encrypt_decrypt() {
        let text = "圣光会抛弃你的，英雄，就像抛弃我那样。——巫妖王";

        let prk = "6aea1ccf610488aaa7fddba3dd6d76d3bdfd50f957d847be3d453defb695f28e";
        let puk = "04a8af64e38eea41c254df769b5b41fbaa2d77b226b301a2636d463c52b46c777230ad1714e686dd641b9e04596530b38f6a64215b0ed3b081f8641724c5443a6e";

        let cipher = sm2::encrypt(puk, text);
        let plain = sm2::decrypt(prk, &cipher);

        assert_eq!(plain, text);
    }

    #[test]
    fn sm2_sign_verify() {
        let text = "圣光会抛弃你的，英雄，就像抛弃我那样。——巫妖王";

        let prk = "6aea1ccf610488aaa7fddba3dd6d76d3bdfd50f957d847be3d453defb695f28e";
        let puk = "04a8af64e38eea41c254df769b5b41fbaa2d77b226b301a2636d463c52b46c777230ad1714e686dd641b9e04596530b38f6a64215b0ed3b081f8641724c5443a6e";

        let s = sm2::sign(prk, puk, text);
        let f = sm2::verify(puk, text, &s);

        assert_eq!(f, true);
    }

    #[test]
    fn sm3_hash() {
        let hash = sm3::digest("abc");
        assert_eq!(hash, "66c7f0f462eeedd9d1f2d46bdc10e4e24167c4875cf2f7a2297da02b8f4ba8e0");
    }


    #[test]
    fn sm4_key() {
        let key = sm4::generate_key();
        println!("key = {:?}", key);
    }

    #[test]
    fn sm4_iv() {
        let iv = sm4::generate_iv();
        println!("iv = {:?}", iv);
    }

    #[test]
    fn sm4_ecb() {
        let key = sm4::generate_key();
        let plain = "圣光会抛弃你的，英雄，就像抛弃我那样。——巫妖王";
        let mode = sm4::Mode::ECB { key };

        let crypto = sm4::CryptoFactory::new(mode);
        // 加密
        let cipher = crypto.encrypt(String::from(plain));
        // 解密
        let text = crypto.decrypt(cipher);

        assert_eq!(plain, text);
    }

    #[test]
    fn sm4_cbc() {
        let key = sm4::generate_key();
        let iv = sm4::generate_iv();
        let plain = "记住‘被遗忘者’的含义，我们既非生者也非死者，我们将被活着的和死去的人遗忘。\
        我们回到了曾经告别的世界上，但是却永远无法回到我们曾经活着的那些日子，永远无法回到那些我们曾经爱过的人的身边。\
        我们是存在也是诅咒，因此我们遗忘过去，并且被过去遗忘。——希尔瓦娜斯";
        let mode = sm4::Mode::CBC { key, iv };

        let crypto = sm4::CryptoFactory::new(mode);
        // 加密
        let cipher = crypto.encrypt(String::from(plain));
        // 解密
        let text = crypto.decrypt(cipher);
        assert_eq!(plain, text);
    }

    #[test]
    fn sm4_cfb() {
        let key = sm4::generate_key();
        let iv = sm4::generate_iv();
        let plain = "兽人永不为奴，我们终将成王。——加尔鲁什·地狱咆哮";
        let mode = sm4::Mode::CFB { key, iv };

        let crypto = sm4::CryptoFactory::new(mode);
        // 加密
        let cipher = crypto.encrypt(String::from(plain));
        // 解密
        let text = crypto.decrypt(cipher);
        assert_eq!(plain, text);
    }

    #[test]
    fn sm4_ofb() {
        let key = sm4::generate_key();
        let iv = sm4::generate_iv();
        let plain = "没人生来杰出（No one breather who is worthier）——奥格瑞姆·毁灭之锤";
        let mode = sm4::Mode::OFB { key, iv };

        let crypto = sm4::CryptoFactory::new(mode);
        // 加密
        let cipher = crypto.encrypt(String::from(plain));
        // 解密
        let text = crypto.decrypt(cipher);
        assert_eq!(plain, text);
    }

    #[test]
    fn sm4_ctr() {
        let key = sm4::generate_key();
        let iv = sm4::generate_iv();
        let plain = "有时候…我真希望能有人来给我一个又大又久的拥抱。  ——克尔苏加德";
        let mode = sm4::Mode::CTR { key, iv };

        let crypto = sm4::CryptoFactory::new(mode);
        // 加密
        let cipher = crypto.encrypt(String::from(plain));
        // 解密
        let text = crypto.decrypt(cipher);
        assert_eq!(plain, text);
    }
}