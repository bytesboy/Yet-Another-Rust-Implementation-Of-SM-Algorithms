// use mlua::Lua;
// use mlua::prelude::{LuaResult, LuaTable};
//
// mod sm2 {
//     use mlua::Lua;
//     use mlua::prelude::LuaResult;
//
//     pub(crate) fn generate_key(_: &Lua, _: ()) -> LuaResult<(String, String)> {
//         let pair = crate::sm2::generate_key();
//         Ok(pair)
//     }
//
//     pub(crate) fn encrypt(_: &Lua, (public_key, plain): (String, String)) -> LuaResult<String> {
//         let cipher = crate::sm2::encrypt(&public_key, &plain);
//         Ok(cipher)
//     }
//
//     pub(crate) fn decrypt(_: &Lua, (private_key, cipher): (String, String)) -> LuaResult<String> {
//         let plain = crate::sm2::decrypt(&private_key, &cipher);
//         Ok(plain)
//     }
//
//     pub(crate) fn encrypt_c1c2c3(_: &Lua, (public_key, plain): (String, String)) -> LuaResult<String> {
//         let cipher = crate::sm2::encrypt_c1c2c3(&public_key, &plain);
//         Ok(cipher)
//     }
//
//     pub(crate) fn decrypt_c1c2c3(_: &Lua, (private_key, cipher): (String, String)) -> LuaResult<String> {
//         let plain = crate::sm2::decrypt_c1c2c3(&private_key, &cipher);
//         Ok(plain)
//     }
//
//     pub(crate) fn sign(_: &Lua, (private_key, public_key, plain): (String, String, String)) -> LuaResult<String> {
//         let signature = crate::sm2::sign(&private_key, &public_key, &plain);
//         Ok(signature)
//     }
//
//     pub(crate) fn verify(_: &Lua, (public_key, plain, signature): (String, String, String)) -> LuaResult<bool> {
//         let flag = crate::sm2::verify(&public_key, &plain, &signature);
//         Ok(flag)
//     }
// }
//
//
// mod sm3 {
//     use mlua::Lua;
//     use mlua::prelude::LuaResult;
//
//     pub(crate) fn digest(_: &Lua, data: String) -> LuaResult<String> {
//         let hash = crate::sm3::digest(&data);
//         Ok(hash)
//     }
// }
//
//
// mod sm4 {
//     use mlua::Lua;
//     use mlua::prelude::LuaResult;
//
//
//     pub(crate) fn generate_key(_: &Lua, _: ()) -> LuaResult<String> {
//         let key = crate::sm4::generate_key();
//         Ok(key)
//     }
//
//     pub(crate) fn generate_iv(_: &Lua, _: ()) -> LuaResult<String> {
//         let iv = crate::sm4::generate_iv();
//         Ok(iv)
//     }
//
//     pub(crate) fn encrypt_ecb(_: &Lua, (key, plain): (String, String)) -> LuaResult<String> {
//         let cipher = crate::sm4::encrypt_ecb(key, plain);
//         Ok(cipher)
//     }
//
//     pub(crate) fn decrypt_ecb(_: &Lua, (key, cipher): (String, String)) -> LuaResult<String> {
//         let plain = crate::sm4::decrypt_ecb(key, cipher);
//         Ok(plain)
//     }
//
//     pub(crate) fn encrypt_cbc(_: &Lua, (key, iv, plain): (String, String, String)) -> LuaResult<String> {
//         let cipher = crate::sm4::encrypt_cbc(key, iv, plain);
//         Ok(cipher)
//     }
//
//     pub(crate) fn decrypt_cbc(_: &Lua, (key, iv, cipher): (String, String, String)) -> LuaResult<String> {
//         let plain = crate::sm4::decrypt_cbc(key, iv, cipher);
//         Ok(plain)
//     }
//
//     pub(crate) fn encrypt_cfb(_: &Lua, (key, iv, plain): (String, String, String)) -> LuaResult<String> {
//         let cipher = crate::sm4::encrypt_cfb(key, iv, plain);
//         Ok(cipher)
//     }
//
//     pub(crate) fn decrypt_cfb(_: &Lua, (key, iv, cipher): (String, String, String)) -> LuaResult<String> {
//         let plain = crate::sm4::decrypt_cfb(key, iv, cipher);
//         Ok(plain)
//     }
//
//     pub(crate) fn encrypt_ofb(_: &Lua, (key, iv, plain): (String, String, String)) -> LuaResult<String> {
//         let cipher = crate::sm4::encrypt_ofb(key, iv, plain);
//         Ok(cipher)
//     }
//
//     pub(crate) fn decrypt_ofb(_: &Lua, (key, iv, cipher): (String, String, String)) -> LuaResult<String> {
//         let plain = crate::sm4::decrypt_ofb(key, iv, cipher);
//         Ok(plain)
//     }
//
//     pub(crate) fn encrypt_ctr(_: &Lua, (key, iv, plain): (String, String, String)) -> LuaResult<String> {
//         let cipher = crate::sm4::encrypt_ctr(key, iv, plain);
//         Ok(cipher)
//     }
//
//     pub(crate) fn decrypt_ctr(_: &Lua, (key, iv, cipher): (String, String, String)) -> LuaResult<String> {
//         let plain = crate::sm4::decrypt_ctr(key, iv, cipher);
//         Ok(plain)
//     }
// }
//
//
// #[mlua::lua_module]
// fn libyarism(lua: &Lua) -> LuaResult<LuaTable> {
//     let exports = lua.create_table()?;
//
//     exports.set("sm2_generate_key", lua.create_function(sm2::generate_key)?)?;
//     exports.set("sm2_encrypt", lua.create_function(sm2::encrypt)?)?;
//     exports.set("sm2_decrypt", lua.create_function(sm2::decrypt)?)?;
//     exports.set("sm2_encrypt_c1c2c3", lua.create_function(sm2::encrypt_c1c2c3)?)?;
//     exports.set("sm2_decrypt_c1c2c3", lua.create_function(sm2::decrypt_c1c2c3)?)?;
//     exports.set("sm2_sign", lua.create_function(sm2::sign)?)?;
//     exports.set("sm2_verify", lua.create_function(sm2::verify)?)?;
//
//     exports.set("sm3_digest", lua.create_function(sm3::digest)?)?;
//
//     exports.set("sm4_generate_key", lua.create_function(sm4::generate_key)?)?;
//     exports.set("sm4_generate_iv", lua.create_function(sm4::generate_iv)?)?;
//     exports.set("sm4_encrypt_ecb", lua.create_function(sm4::encrypt_ecb)?)?;
//     exports.set("sm4_decrypt_ecb", lua.create_function(sm4::decrypt_ecb)?)?;
//     exports.set("sm4_encrypt_cbc", lua.create_function(sm4::encrypt_cbc)?)?;
//     exports.set("sm4_decrypt_cbc", lua.create_function(sm4::decrypt_cbc)?)?;
//     exports.set("sm4_encrypt_cfb", lua.create_function(sm4::encrypt_cfb)?)?;
//     exports.set("sm4_decrypt_cfb", lua.create_function(sm4::decrypt_cfb)?)?;
//     exports.set("sm4_encrypt_ofb", lua.create_function(sm4::encrypt_ofb)?)?;
//     exports.set("sm4_decrypt_ofb", lua.create_function(sm4::decrypt_ofb)?)?;
//     exports.set("sm4_encrypt_ctr", lua.create_function(sm4::encrypt_ctr)?)?;
//     exports.set("sm4_decrypt_ctr", lua.create_function(sm4::decrypt_ctr)?)?;
//
//     Ok(exports)
// }
