use std::cmp::Ordering;
use std::mem;
use std::ops::{Add, Mul, Shl, Shr};
use std::sync::Once;

use num_bigint::{BigInt, BigUint, ToBigInt};
use num_integer::Integer;
use num_traits::FromPrimitive;

use crate::sm2::core::{Elliptic, EllipticProvider};

const EC_P: [u8; 32] = [
    0xFF, 0xFF, 0xFF, 0xFE,
    0xFF, 0xFF, 0xFF, 0xFF,
    0xFF, 0xFF, 0xFF, 0xFF,
    0xFF, 0xFF, 0xFF, 0xFF,
    0xFF, 0xFF, 0xFF, 0xFF,
    0x00, 0x00, 0x00, 0x00,
    0xFF, 0xFF, 0xFF, 0xFF,
    0xFF, 0xFF, 0xFF, 0xFF,
];

const EC_A: [u8; 32] = [
    0xFF, 0xFF, 0xFF, 0xFE,
    0xFF, 0xFF, 0xFF, 0xFF,
    0xFF, 0xFF, 0xFF, 0xFF,
    0xFF, 0xFF, 0xFF, 0xFF,
    0xFF, 0xFF, 0xFF, 0xFF,
    0x00, 0x00, 0x00, 0x00,
    0xFF, 0xFF, 0xFF, 0xFF,
    0xFF, 0xFF, 0xFF, 0xFC,
];

const EC_B: [u8; 32] = [
    0x28, 0xE9, 0xFA, 0x9E,
    0x9D, 0x9F, 0x5E, 0x34,
    0x4D, 0x5A, 0x9E, 0x4B,
    0xCF, 0x65, 0x09, 0xA7,
    0xF3, 0x97, 0x89, 0xF5,
    0x15, 0xAB, 0x8F, 0x92,
    0xDD, 0xBC, 0xBD, 0x41,
    0x4D, 0x94, 0x0E, 0x93,
];

const EC_N: [u8; 32] = [
    0xFF, 0xFF, 0xFF, 0xFE,
    0xFF, 0xFF, 0xFF, 0xFF,
    0xFF, 0xFF, 0xFF, 0xFF,
    0xFF, 0xFF, 0xFF, 0xFF,
    0x72, 0x03, 0xDF, 0x6B,
    0x21, 0xC6, 0x05, 0x2B,
    0x53, 0xBB, 0xF4, 0x09,
    0x39, 0xD5, 0x41, 0x23,
];

const EC_GX: [u8; 32] = [
    0x32, 0xC4, 0xAE, 0x2C,
    0x1F, 0x19, 0x81, 0x19,
    0x5F, 0x99, 0x04, 0x46,
    0x6A, 0x39, 0xC9, 0x94,
    0x8F, 0xE3, 0x0B, 0xBF,
    0xF2, 0x66, 0x0B, 0xE1,
    0x71, 0x5A, 0x45, 0x89,
    0x33, 0x4C, 0x74, 0xC7,
];

const EC_GY: [u8; 32] = [
    0xBC, 0x37, 0x36, 0xA2,
    0xF4, 0xF6, 0x77, 0x9C,
    0x59, 0xBD, 0xCE, 0xE3,
    0x6B, 0x69, 0x21, 0x53,
    0xD0, 0xA9, 0x87, 0x7C,
    0xC6, 0x2A, 0x47, 0x40,
    0x02, 0xDF, 0x32, 0xE5,
    0x21, 0x39, 0xF0, 0xA0,
];


/// Montgomery precomputed (1/R) mod P
/// R = 2^257 = 2 * 16^64 = 0x20000000000000000000000000000000000000000000000000000000000000000
const RI: [u8; 32] = [
    0x7F, 0xFF, 0xFF, 0xFD,
    0x80, 0x00, 0x00, 0x02,
    0xFF, 0xFF, 0xFF, 0xFE,
    0x00, 0x00, 0x00, 0x01,
    0x7F, 0xFF, 0xFF, 0xFE,
    0x80, 0x00, 0x00, 0x03,
    0x7F, 0xFF, 0xFF, 0xFC,
    0x80, 0x00, 0x00, 0x02
];

/// MULTI_BASE_POINT_PRECOMPUTED contains precomputed values to aid the calculation of scalar
/// multiples of the base point, G. It's actually two, equal length, tables concatenated.
///
/// The first table contains (x,y) field element pairs for 16 multiples of the base point, G.
///
///   Index  |  Index (binary) | Value
///
///       0  |           0000  | 0G (all zeros, omitted)
///       1  |           0001  | G
///       2  |           0010  | 2**64G
///       3  |           0011  | 2**64G + G
///       4  |           0100  | 2**128G
///       5  |           0101  | 2**128G + G
///       6  |           0110  | 2**128G + 2**64G
///       7  |           0111  | 2**128G + 2**64G + G
///       8  |           1000  | 2**192G
///       9  |           1001  | 2**192G + G
///      10  |           1010  | 2**192G + 2**64G
///      11  |           1011  | 2**192G + 2**64G + G
///      12  |           1100  | 2**192G + 2**128G
///      13  |           1101  | 2**192G + 2**128G + G
///      14  |           1110  | 2**192G + 2**128G + 2**64G
///      15  |           1111  | 2**192G + 2**128G + 2**64G + G
///
/// The second table follows the same style, but the terms are 2^32G, 2^96G, 2^160G, 2^224G.
///
///      16  |          10000  | 2**32G
///      17  |          10010  | 2**96G
///      18  |          10001  | 2**96G  + 2**32G
///      19  |          10011  | 2**160G
///      20  |          10100  | 2**160G + 2**32G
///      21  |          10101  | 2**160G + 2**96G
///      22  |          10110  | 2**160G + 2**96G + 2**32G
///      23  |          10111  | 2**224G
///      24  |          11000  | 2**224G + 2**32G
///      25  |          11001  | 2**224G + 2**96G
///      26  |          11011  | 2**224G + 2**96G + 2**32G
///      27  |          11100  | 2**224G + 2**160G
///      28  |          11101  | 2**224G + 2**160G  + 2**32G
///      29  |          11110  | 2**224G + 2**160G + 2**96G
///      30  |          11111  | 2**224G + 2**160G + 2**96G + 2**32G
///
/// precompute(1) => \[u32; 15 * 9 * 2]
/// precompute(2**32) => \[u32; 15 * 9 * 2]
/// MULTI_BASE_POINT_PRECOMPUTED = \[precompute(1), precompute(2**32)]
const MULTI_BASE_POINT_PRECOMPUTED: [u32; 15 * 2 * 9 * 2] = [
    0x0830053D, 0x0328990F, 0x06C04FE1, 0x0C0F72E5, 0x01E19F3C, 0x0666B093, 0x0175A87B, 0x0EC38276, 0x0222CF4B,
    0x185A1BBA, 0x0354E593, 0x1295FAC1, 0x0F2BC469, 0x047C60FA, 0x0C19B8A9, 0x0F63533E, 0x0903AE6B, 0x0C79ACBA,
    0x15B061A4, 0x033E020B, 0x0DFFB34B, 0x00FCF2C8, 0x16582E08, 0x0262F203, 0x0FB34381, 0x00A55452, 0x0604F0FF,
    0x041F1F90, 0x0D64CED2, 0x0EE377BF, 0x075F05F0, 0x189467AE, 0x00E2244E, 0x1E7700E8, 0x03FBC464, 0x09612D2E,
    0x1341B3B8, 0x0EE84E23, 0x1EDFA5B4, 0x014E6030, 0x19E87BE9, 0x092F533C, 0x1665D96C, 0x0226653E, 0x0A238D3E,
    0x00F5C62C, 0x0095BB7A, 0x1F0E5A41, 0x028789C3, 0x1F251D23, 0x08726609, 0x0E918910, 0x08096848, 0x0F63D028,
    0x152296A1, 0x09F561A8, 0x14D376FB, 0x0898788A, 0x061A95FB, 0x0A59466D, 0x159A003D, 0x01AD1698, 0x093CCA08,
    0x1B314662, 0x0706E006, 0x11CE1E30, 0x0097B710, 0x172FBC0D, 0x08F50158, 0x11C7FFE7, 0x0D182CCE, 0x0C6AD9E8,
    0x12EA31B2, 0x0C4E4F38, 0x175B0D96, 0x0EC06337, 0x075A9C12, 0x0B001FDF, 0x093E82F5, 0x034607DE, 0x0B8035ED,
    0x17F97924, 0x075CF9E6, 0x0DCEAEDD, 0x02529924, 0x1A10C5FF, 0x0B1A54DC, 0x019464D8, 0x002D1997, 0x0DE6A110,
    0x1E276EE5, 0x095C510C, 0x1ACA7C7A, 0x0FE48ACA, 0x121AD4D9, 0x0E4132C6, 0x08239B9D, 0x040EA9CD, 0x00816C7B,
    0x0632D7A4, 0x0A679813, 0x05911FCF, 0x082B0F7C, 0x057B0AD5, 0x000BEF65, 0x0D541365, 0x07F9921F, 0x00C62E7A,
    0x03F4B32D, 0x058E50E1, 0x06427AED, 0x0DCDDA67, 0x0E8C2D3E, 0x06AA54A4, 0x18DF4C35, 0x049A6A8E, 0x03CD3D0C,
    0x00D7ADF2, 0x00CBCA97, 0x1BDA5F2D, 0x03258579, 0x0606B1E6, 0x06FC1B5B, 0x1AC27317, 0x0503CA16, 0x0A677435,
    0x0057BC73, 0x03992A42, 0x0BAB987B, 0x0FAB25EB, 0x128912A4, 0x090A1DC4, 0x1402D591, 0x09FFBCFC, 0x0AA48856,
    0x07A7C2DC, 0x0CEFD08A, 0x1B29BDA6, 0x0A785641, 0x16462D8C, 0x076241B7, 0x079B6C3B, 0x0204AE18, 0x0F41212B,
    0x1F567A4D, 0x0D6CE6DB, 0x0EDF1784, 0x0111DF34, 0x085D7955, 0x055FC189, 0x1B7AE265, 0x0F9281AC, 0x0DED7740,
    0x0F19468B, 0x083763BB, 0x08FF7234, 0x03DA7DF8, 0x09590AC3, 0x0DC96F2A, 0x16E44896, 0x07931009, 0x099D5ACC,
    0x10F7B842, 0x0AEF5E84, 0x0C0310D7, 0x0DEBAC2C, 0x02A7B137, 0x04342344, 0x19633649, 0x03A10624, 0x04B4CB56,
    0x1D809C59, 0x00AC007F, 0x1F0F4BCD, 0x0A1AB06E, 0x0C5042CF, 0x082C0C77, 0x076C7563, 0x022C30F3, 0x03BF1568,
    0x07A895BE, 0x0FCCA554, 0x12E90E4C, 0x07B4AB5F, 0x13AEB76B, 0x05887E2C, 0x1D7FE1E3, 0x0908C8E3, 0x095800EE,
    0x0B36BD54, 0x0F08905D, 0x04E73AE8, 0x0F5A7E48, 0x00A67CB0, 0x050E1067, 0x1B944A0A, 0x0F29C83A, 0x0B23CFB9,
    0x07A895BE, 0x0FCCA554, 0x12E90E4C, 0x07B4AB5F, 0x13AEB76B, 0x05887E2C, 0x1D7FE1E3, 0x0908C8E3, 0x095800EE,
    0x0B36BD54, 0x0F08905D, 0x04E73AE8, 0x0F5A7E48, 0x00A67CB0, 0x050E1067, 0x1B944A0A, 0x0F29C83A, 0x0B23CFB9,
    0x11237F01, 0x0E2A820B, 0x0FD53B95, 0x06BEB5EE, 0x1AEB790C, 0x0E470D53, 0x02C2CFEE, 0x01C1D8D8, 0x0A520FC4,
    0x1518E034, 0x0A584DD4, 0x029E572B, 0x0D4594FC, 0x141A8F6F, 0x08DFCCF3, 0x05D20BA3, 0x02EB60C3, 0x09F16EB0,
    0x11CEC356, 0x0F039F84, 0x1B0990C1, 0x0C91E526, 0x10B65BAE, 0x0F0616E8, 0x173FA3FF, 0x0EC8CCF9, 0x0BE32790,
    0x11DA3E79, 0x0E2F35C7, 0x0908875C, 0x0DACF7BD, 0x0538C165, 0x08D1487F, 0x07C31AED, 0x021AF228, 0x07E1689D,
    0x0DFC23CA, 0x024F15DC, 0x025EF3C4, 0x035248CD, 0x099A0F43, 0x0A4B6ECC, 0x00D066B3, 0x02481152, 0x037A7688,
    0x15A444B6, 0x0B62300C, 0x004B841B, 0x0A655E79, 0x0D53226D, 0x0BEB348A, 0x0127F3C2, 0x0B989247, 0x071A277D,
    0x19E9DFCB, 0x0B8F92D0, 0x0E2D226C, 0x0390A8B0, 0x183CC462, 0x07BD8167, 0x1F32A552, 0x05E02DB4, 0x0A146EE9,
    0x1A003957, 0x01C95F61, 0x1EEEC155, 0x026F811F, 0x0F9596BA, 0x03082BFB, 0x096DF083, 0x03E3A289, 0x07E2D8BE,
    0x157A63E0, 0x099B8941, 0x1DA7D345, 0x00CC6CD0, 0x10BEED9A, 0x048E83C0, 0x13AA2E25, 0x07CAD710, 0x04029988,
    0x13DFA9DD, 0x0B94F884, 0x1F4ADFEF, 0x00B88543, 0x16F5F8DC, 0x0A6A67F4, 0x14E274E2, 0x05E56CF4, 0x002F24EF,
    0x1E9EF967, 0x0FE09BAD, 0x0FE079B3, 0x0CC0AE9E, 0x0B3EDF6D, 0x03E961BC, 0x130D7831, 0x031043D6, 0x0BA986F9,
    0x01D28055, 0x065240CA, 0x04971FA3, 0x081B17F8, 0x11EC34A5, 0x08366DDC, 0x01471809, 0x0FA5F1C6, 0x0C911E15,
    0x08849491, 0x0CF4C2E2, 0x14471B91, 0x039F75BE, 0x0445C21E, 0x0F1585E9, 0x072CC11F, 0x04C79F0C, 0x0E5522E1,
    0x1874C1EE, 0x04444211, 0x07914884, 0x03D1B133, 0x0025BA3C, 0x04194F65, 0x1C0457EF, 0x0AC4899D, 0x0E1FA66C,
    0x130A7918, 0x09B8D312, 0x04B1C5C8, 0x061CCAC3, 0x18C8AA6F, 0x0E93CB0A, 0x0DCCB12C, 0x0DE10825, 0x0969737D,
    0x0F58C0C3, 0x07CEE6A9, 0x0C2C329A, 0x0C7F9ED9, 0x107B3981, 0x0696A40E, 0x152847FF, 0x04D88754, 0x0B141F47,
    0x05A16FFE, 0x03A7870A, 0x18667659, 0x03B72B03, 0x0B1C9435, 0x09285394, 0x0A00005A, 0x0037506C, 0x02EDC0BB,
    0x19AFE392, 0x0EB39CAC, 0x177EF286, 0x0DF87197, 0x19F844ED, 0x00031FE8, 0x015F9BFD, 0x0080DBEC, 0x0342E96E,
    0x0497ACED, 0x0E88E909, 0x1F5FA9BA, 0x0530A6EE, 0x1EF4E3F1, 0x069FFD12, 0x0583006D, 0x02ECC9B1, 0x0362DB70,
    0x18C7BDC5, 0x0F4BB3C5, 0x1C90B957, 0x0F067C09, 0x09768F2B, 0x0F73566A, 0x1939A900, 0x0198C38A, 0x0202A2A1,
    0x04BBF5A6, 0x04E265BC, 0x1F44B6E7, 0x0185CA49, 0x0A39E81B, 0x024AFF5B, 0x04ACC9C2, 0x0638BDD3, 0x0B65B2A8,
    0x06DEF8BE, 0x0B94537A, 0x10B81DEE, 0x0E00EC55, 0x02F2CDF7, 0x0C20622D, 0x02D20F36, 0x0E03C8C9, 0x0898EA76,
    0x08E3921B, 0x08905BFF, 0x1E94B6C8, 0x0EE7AD86, 0x154797F2, 0x0A620863, 0x03FBD0D9, 0x01F3CAAB, 0x030C24BD,
    0x19D3892F, 0x059C17A2, 0x1AB4B0AE, 0x0F8714EE, 0x090C4098, 0x0A9C800D, 0x1910236B, 0x0EA808D3, 0x09AE2F31,
    0x1A15AD64, 0x0A48C8D1, 0x184635A4, 0x0B725EF1, 0x11921DCC, 0x03F866DF, 0x16C27568, 0x0BDF580A, 0x0B08F55C,
    0x0186EE1C, 0x0B1627FA, 0x034E82F6, 0x0933837E, 0x0F311BE5, 0x0FEDB03B, 0x167F72CD, 0x0A5469C0, 0x09C82531,
    0x0B92A24B, 0x014FDC8B, 0x141980D1, 0x0BDC3A49, 0x07E02BB1, 0x0AF4E6DD, 0x106D99E1, 0x0D4616FC, 0x093C2717,
    0x1C0A0507, 0x0C6D5FED, 0x09A03D8B, 0x0A1D22B0, 0x127853E3, 0x0C4AC6B8, 0x1A048CF7, 0x09AFB72C, 0x065D485D,
    0x0B92A24B, 0x014FDC8B, 0x141980D1, 0x0BDC3A49, 0x07E02BB1, 0x0AF4E6DD, 0x106D99E1, 0x0D4616FC, 0x093C2717,
    0x1C0A0507, 0x0C6D5FED, 0x09A03D8B, 0x0A1D22B0, 0x127853E3, 0x0C4AC6B8, 0x1A048CF7, 0x09AFB72C, 0x065D485D,
    0x0175BCBB, 0x0B29B49F, 0x1806B79C, 0x012FB61F, 0x170B3A10, 0x03AAF1CF, 0x0A224085, 0x079D26AF, 0x097759E2,
    0x092E19F1, 0x0B32714D, 0x1F00D9F1, 0x0C728619, 0x09E6F627, 0x0E745E24, 0x18EA4ACE, 0x0FC60A41, 0x0125F5B2,
    0x0C3CF512, 0x039ED486, 0x0F4D15FA, 0x0F9167FD, 0x1C1F5DD5, 0x0C21A53E, 0x01897930, 0x0957A112, 0x021059A0,
    0x1F9E3DDC, 0x0A4DFCED, 0x08427F6F, 0x0726FBE7, 0x1EA658F8, 0x02FDCD4C, 0x17E9B66F, 0x0B2E7C2E, 0x039923BF,
    0x01BAE104, 0x03973CE5, 0x0C6F264C, 0x03511B84, 0x124195D7, 0x011996BD, 0x020BE23D, 0x0DC437C4, 0x04B4F16B,
    0x011902A0, 0x06C29CC9, 0x1D5FFBE6, 0x0DB0B4C7, 0x10144C14, 0x02F2B719, 0x00301189, 0x02343336, 0x0A0BF2AC,
];

const P256CARRY: [u32; 8 * 9] = [
    0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
    0x00000002, 0x00000000, 0x1FFFFF00, 0x000007FF, 0x00000000, 0x00000000, 0x00000000, 0x02000000, 0x00000000,
    0x00000004, 0x00000000, 0x1FFFFE00, 0x00000FFF, 0x00000000, 0x00000000, 0x00000000, 0x04000000, 0x00000000,
    0x00000006, 0x00000000, 0x1FFFFD00, 0x000017FF, 0x00000000, 0x00000000, 0x00000000, 0x06000000, 0x00000000,
    0x00000008, 0x00000000, 0x1FFFFC00, 0x00001FFF, 0x00000000, 0x00000000, 0x00000000, 0x08000000, 0x00000000,
    0x0000000A, 0x00000000, 0x1FFFFB00, 0x000027FF, 0x00000000, 0x00000000, 0x00000000, 0x0A000000, 0x00000000,
    0x0000000C, 0x00000000, 0x1FFFFA00, 0x00002FFF, 0x00000000, 0x00000000, 0x00000000, 0x0C000000, 0x00000000,
    0x0000000E, 0x00000000, 0x1FFFF900, 0x000037FF, 0x00000000, 0x00000000, 0x00000000, 0x0E000000, 0x00000000,
];

const P256ZERO31: [u32; 9] = [
    0x7FFFFFF8, 0x3FFFFFFC, 0x800003FC, 0x3FFFDFFC, 0x7FFFFFFC, 0x3FFFFFFC, 0x7FFFFFFC, 0x37FFFFFC, 0x7FFFFFFC
];

const P256FACTOR: [[u32; 9]; 9] = [
    [0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000],
    [0x00000002, 0x00000000, 0x1FFFFF00, 0x000007FF, 0x00000000, 0x00000000, 0x00000000, 0x02000000, 0x00000000],
    [0x00000004, 0x00000000, 0x1FFFFE00, 0x00000FFF, 0x00000000, 0x00000000, 0x00000000, 0x04000000, 0x00000000],
    [0x00000006, 0x00000000, 0x1FFFFD00, 0x000017FF, 0x00000000, 0x00000000, 0x00000000, 0x06000000, 0x00000000],
    [0x00000008, 0x00000000, 0x1FFFFC00, 0x00001FFF, 0x00000000, 0x00000000, 0x00000000, 0x08000000, 0x00000000],
    [0x0000000A, 0x00000000, 0x1FFFFB00, 0x000027FF, 0x00000000, 0x00000000, 0x00000000, 0x0A000000, 0x00000000],
    [0x0000000C, 0x00000000, 0x1FFFFA00, 0x00002FFF, 0x00000000, 0x00000000, 0x00000000, 0x0C000000, 0x00000000],
    [0x0000000E, 0x00000000, 0x1FFFF900, 0x000037FF, 0x00000000, 0x00000000, 0x00000000, 0x0E000000, 0x00000000],
    [0x00000010, 0x00000000, 0x1FFFF800, 0x00003FFF, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000001],
];


#[derive(Clone, Debug)]
pub struct P256Elliptic {
    ec: Elliptic,
    ri: BigUint,
}

impl P256Elliptic {
    pub fn init() -> Self {
        static mut ELLIPTIC: *const P256Elliptic = std::ptr::null::<P256Elliptic>();
        static INITIALIZER: Once = Once::new();
        unsafe {
            INITIALIZER.call_once(|| {
                let p256 = P256Elliptic {
                    ec: Elliptic {
                        p: BigUint::from_bytes_be(&EC_P),
                        a: BigUint::from_bytes_be(&EC_A),
                        b: BigUint::from_bytes_be(&EC_B),
                        gx: BigUint::from_bytes_be(&EC_GX),
                        gy: BigUint::from_bytes_be(&EC_GY),
                        n: BigUint::from_bytes_be(&EC_N),
                        bits: 256,
                    },
                    ri: BigUint::from_bytes_be(&RI),
                };
                ELLIPTIC = mem::transmute(Box::new(p256));
            });
            (*ELLIPTIC).clone()
        }
    }
}

impl EllipticProvider for P256Elliptic {
    fn blueprint(&self) -> &Elliptic {
        &self.ec
    }

    fn scalar_multiply(&self, x: BigUint, y: BigUint, k: BigUint) -> (BigUint, BigUint) {
        let point = P256Point(x, y);
        let p = point.multiply(k);
        (p.0, p.1)
    }

    fn scalar_base_multiply(&self, k: BigUint) -> (BigUint, BigUint) {
        let elliptic = self.ec.clone();
        let base = P256BasePoint {
            point: P256Point(elliptic.gx.clone(), elliptic.gy.clone()),
            order: elliptic.n,
        };
        let p = base.multiply(k);
        (p.0, p.1)
    }
}

trait Multiplication {
    fn multiply(&self, scalar: BigUint) -> P256Point;
}

/// Jacobian coordinates: (x, y, z)  y^2 = x^3 + axz^4 + bz^6
/// Affine coordinates: (X = x/z^2, Y = y/z^3)  Y^2 = X^3 + aX +b
#[derive(Clone, Debug)]
struct P256Point(BigUint, BigUint);

/// 基点
#[derive(Clone, Debug)]
struct P256BasePoint {
    point: P256Point,
    order: BigUint,
}

impl Multiplication for P256Point {
    fn multiply(&self, scalar: BigUint) -> P256Point {
        todo!()
    }
}

impl Multiplication for P256BasePoint {
    fn multiply(&self, scalar: BigUint) -> P256Point {
        let scalar = {
            // compare scalar and order, n = (scalar mod order) if scalar > order else scalar
            if let Ordering::Greater = scalar.cmp(&self.order) {
                scalar.mod_floor(&self.order)
            } else {
                scalar
            }
        };

        let mut scalar_bytes = [0u8; 32];
        for (i, v) in scalar.to_bytes_le().iter().enumerate() {
            scalar_bytes[i] = *v;
        }


        self.point.clone()
    }
}


/// Field elements are represented as nine, unsigned 32-bit words. The value of a field element is:
///
/// ```
/// Value = (x8 * 2^228) + (x7 * 2^200) + (x6 * 2^171) + (x5 * 2^143) + (x4 * 2^114) + (x3 * 2^86) +
///         (x2 * 2^57)  + (x1 * 2^29)  + x0
/// ```
///
/// That is, each limb is alternately 29 or 28-bits wide in little-endian order.
///
/// This means that a field element hits 2^257, rather than 2^256 as we would like.
/// A 28, 29, ... pattern would cause us to hit 2^256, but that causes problems
/// when multiplying as terms end up one bit short of a limb
/// which would require much bit-shifting to correct.
///
/// ***
/// ### Pattern
/// * bits-257: |29bits|28bits|29bits|29bits|28bits|29bits|29bits|28bits|29bits|
///
/// Finally, the values stored in a field element are in Montgomery form.
/// So the value |y| is stored as (y*R) mod p, where p is the P-256 prime and R is 2^257.
/// ***
/// ### little-endian order
///
/// Example: payload = \[x0, x1, x2, x3, x4, x5, x6, x7, x8]
///
/// | 29bits | 28bits | 29bits | 28bits | 29bits | 28bits | 29bits | 28bits | 29bits |
/// | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ |
/// |   x8   |   x7   |   x6   |   x5   |   x4   |   x3   |   x2   |   x1   |   x0   |
///
/// 0xFFFFFFF  = 1111111111111111111111111111
/// 0x1FFFFFFF = 11111111111111111111111111111
enum LimbPattern {
    WIDTH28BITS = 0xFFFFFFF,
    WIDTH29BITS = 0x1FFFFFFF,
}

struct Payload {
    data: [u32; 9],
}

struct PayloadHelper;

impl PayloadHelper {
    /// ### Example
    ///
    /// n: 115792089210356248756420345214020892766250353991924191454421193933289684991996
    ///
    /// * step 1 :
    ///     ```
    ///     x = (n * 2^257) % p;
    ///     x = 115792089048596568753516506446018802244132569949625955944202853485549017104377
    ///       = 1111111111111111111111111111100011111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111001000000000000000000000000000001101111111111111111111111111111111111111111111111111111111111111001
    ///     ```
    /// * step 2: while loop, every time extract 29bits and 28 bits.
    ///     + step 2.1.1: change x: bigint to x: vec\<u32>
    ///         ```
    ///         x = [4294967289, 4294967295, 6, 4294967289, 4294967295, 4294967295, 4294967295, 4294967288]
    ///         ```
    ///     + step 2.1.2: extract 29 bits of x using operator &, `x[0] & 0x1FFFFFFF`
    ///         ```
    ///         11111111111111111111111111111001
    ///         &  11111111111111111111111111111
    ///         =  11111111111111111111111111001
    ///         ```
    ///     + step 2.1.3: right shift 29 bits on purpose to delete the extracted 29 bits.
    ///        ```
    ///         x = 11111111111111111111111111111000111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111110010000000000000000000000000000011011111111111111111111111111111111111
    ///        ```
    ///     + step 2.2.1: some operation like 2.1.1
    ///     + step 2.2.2: extract 28 bits of x using operator &, `x[0] & 0xFFFFFFF`
    ///     + step 2.2.3: right shift 28 bits.
    /// * step 3: get the result
    ///     ```
    ///     data = [536870905, 268435455, 895, 268428288, 536870911, 268435455, 536870911, 150994943, 268435455]
    fn transform(n: &BigInt) -> Payload {
        let elliptic = P256Elliptic::init();
        let mut data: [u32; 9] = [0; 9];
        // data = n * R mod p = n * 2^257 % p
        let mut x: BigInt = BigInt::shl(n.clone(), 257);
        x = x.mod_floor(&elliptic.ec.p.to_bigint().unwrap());
        let mut i: usize = 0;
        while i < 9 {
            // x -> [u32]
            let bits = x.to_u32_digits().1;
            if !bits.is_empty() {
                // extract 29 bits using operator &
                data[i] = bits[0] & (LimbPattern::WIDTH29BITS as u32);
            } else {
                data[i] = 0
            }
            // right shift 29 bits
            x = BigInt::shr(x, 29);
            i += 1;
            if i == 9 {
                break;
            }
            let bits = x.to_u32_digits().1;
            if !bits.is_empty() {
                // extract 28 bits using operator &
                data[i] = bits[0] & (LimbPattern::WIDTH28BITS as u32);
            } else {
                data[i] = 0
            }
            // right shift 28 bits
            x = BigInt::shr(x, 28);
            i += 1;
        }
        Payload { data }
    }

    /// Example: payload = \[x0, x1, x2, x3, x4, x5, x6, x7, x8]
    ///
    ///           n = x8
    /// * i=7  => n = n * 2^28 + x7 = x8 * 2^28 + x7
    /// * i=6  => n = n * 2^29 + x6 = x8 * 2^57 + x7 * 2^29 + x6
    /// * i=5  => n = n * 2^28 + x5 = x8 * 2^85 + x7 * 2^57 + x6 * 2^28 + x5
    /// * i=4  => n = n * 2^29 + x4 = x8 * 2^114 + x7 * 2^86 + x6 * 2^57 + x5 * 2^29 + x4
    /// * i=3  => n = n * 2^28 + x3 = x8 * 2^142 + x7 * 2^114 + x6 * 2^85 + x5 * 2^57 + x4 * 2^28 + x3
    /// * i=2  => n = n * 2^29 + x2 = x8 * 2^171 + x7 * 2^143 + x6 * 2^114 + x5 * 2^86 + x4 * 2^57 + x3 * 2^29 + x2
    /// * i=1  => n = n * 2^28 + x1 = x8 * 2^199 + x7 * 2^171 + x6 * 2^142 + x5 * 2^114 + x4 * 2^85 + x3 * 2^57 + x2 * 2^28 + x1
    /// * i=0  => n = n * 2^29 + x0 = x8 * 2^228 + x7 * 2^200 + x6 * 2^171 + x5 * 2^143 + x4 * 2^114 + x3 * 2^86 + x2 * 2^57 + x1 * 2^29 + x0
    ///
    /// return n * RI mod p
    fn restore(payload: &Payload) -> BigInt {
        let elliptic = P256Elliptic::init();
        let mut n = BigInt::from_u32(payload.data[8]).unwrap();
        let mut temp: BigInt;
        let mut i: isize = 7;
        while i >= 0 {
            // i & 1 = 0 => i is even, else i is odd
            if (i & 1) == 0 {
                // even index, n * 2^29
                n = n.shl(29);
            } else {
                // odd index, n * 2^28
                n = n.shl(28);
            }
            temp = BigInt::from_u32(payload.data[i as usize]).unwrap();
            n = n.add(temp);
            i -= 1;
        }
        // formula: data = n * R mod P  => n = data * RI mod p
        n = n.mul(elliptic.ri.to_bigint().unwrap());
        n = n.mod_floor(&elliptic.ec.p.to_bigint().unwrap());
        n
    }


    /// reduce_carry adds a multiple of p in order to cancel |carry|,which is a term at 2^257.
    ///
    /// payload = \[r0, r1, r2, r3, r4, r5, r6, r7, r8]
    ///
    /// we can count Res = carry * 2^257 + r8 * 2^228 + r7 * 2^200 + r6 * 2^171 + r5 * 2^143 + r4 * 2^114 + r3 * 2^86 + r2 * 2^57 + r1 * 2^29 + r0,
    /// and carry * 2^257 could transform to array P256CARRY, through with PayloadHelper::transform(carry), where carry is from 0 to 7.
    ///
    /// So we can mark ResArray = payload + P256CARRY, and Res = PayloadHelper::restore(ResArray)
    ///
    /// |capacity| 32bits |    32bits   | 32bits | 32bits | 32bits |     32bits   |    32bits   | 32bits |    32bits   |
    /// | length | 29bits |   <=29bits  | 29bits | 28bits | 29bits |   <=29bits   |   <=30bits  | 28bits |   <=30bits  |
    /// | ------ | ------ | ----------- | ------ | ------ | ------ | -----------  | ----------- | ------ | ----------- |
    /// |        |   r8   | r7+T[c*9+7] |   r6   |   r5   |   r4   |  r3+T[c*9+3] | r2+T[c*9+2] |   r1   | r0+T[c*9+2] |
    ///
    /// On entry: carry < 2^3, payload\[0,2,...] < 2^29, payload\[1,3,...] < 2^28.
    /// On exit: payload\[0,2,..] < 2^30, payload\[1,3,...] < 2^29.
    fn reduce_carry(payload: &mut Payload, carry: usize) {
        payload.data[0] += P256CARRY[carry * 9 + 0];
        payload.data[2] += P256CARRY[carry * 9 + 2];
        payload.data[3] += P256CARRY[carry * 9 + 3];
        payload.data[7] += P256CARRY[carry * 9 + 7];
    }

    /// reduce_degree sets a = b/R mod p where b contains 64-bit words with the same
    /// 29,28,... bit positions as a field element.
    ///
    /// The values in field elements are in Montgomery form: x*R mod p where R = 2^257.
    /// Since we just multiplied two Montgomery values together, the result is x * y * R * R mod p.
    /// We wish to divide by R in order for the result also to be in Montgomery form.
    ///
    /// On entry: tmp\[i] < 2^64
    /// On exit:  a\[0,2,...] < 2^30, a\[1,3,...] < 2^29
    ///
    /// Limb number:   0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10...
    /// Width (bits):  29| 28| 29| 28| 29| 28| 29| 28| 29| 28| 29
    /// Start bit:     0 | 29| 57| 86|114|143|171|200|228|257|285
    /// (odd phase):   0 | 28| 57| 85|114|142|171|199|228|256|285
    fn reduce_degree(a: &mut Payload, b: &mut [u64; 17]) {
        let mut tmp: [u32; 18] = [0; 18];
        let mut carry: u32;
        let mut x: u32;
        let mut x_mask: u32;

        tmp[0] = (b[0] as u32) & LimbPattern::WIDTH29BITS;
        tmp[1] = (b[0] as u32) >> 29;
        tmp[1] |= (((b[0] >> 32) as u32) << 3) & LimbPattern::WIDTH28BITS;
        tmp[1] += (b[1] as u32) & LimbPattern::WIDTH28BITS;
        carry = tmp[1] >> 28;
        tmp[1] &= LimbPattern::WIDTH28BITS;

        let mut i = 2;
        while i < 17 {
            tmp[i] = ((b[i - 2] >> 32) as u32) >> 25;
            tmp[i] += ((b[i - 1]) as u32) >> 28;
            tmp[i] += (((b[i - 1] >> 32) as u32) << 4) & LimbPattern::WIDTH29BITS;
            tmp[i] += (b[i] as u32) & LimbPattern::WIDTH29BITS;
            tmp[i] += carry;
            carry = tmp[i] >> 29;
            tmp[i] &= LimbPattern::WIDTH29BITS;

            i += 1;
            if i == 17 {
                break;
            }

            tmp[i] = ((b[i - 2] >> 32) as u32) >> 25;
            tmp[i] += (b[i - 1] as u32) >> 29;
            tmp[i] += (((b[i - 1] >> 32) as u32) << 3) & LimbPattern::WIDTH28BITS;
            tmp[i] += (b[i] as u32) & LimbPattern::WIDTH28BITS;
            tmp[i] += carry;
            carry = tmp[i] >> 28;
            tmp[i] &= LimbPattern::WIDTH28BITS;

            i += 1
        }

        tmp[17] = ((b[15] >> 32) as u32) >> 25;
        tmp[17] += (b[16] as u32) >> 29;
        tmp[17] += ((b[16] >> 32) as u32) << 3;
        tmp[17] += carry;

        i = 0;
        loop {
            tmp[i + 1] += tmp[i] >> 29;
            x = tmp[i] & LimbPattern::WIDTH29BITS;
            tmp[i] = 0;

            if x > 0 {
                let mut set4: u32 = 0;
                let mut set7: u32 = 0;
                x_mask = Self::mask(x);
                tmp[i + 2] += (x << 7) & LimbPattern::WIDTH29BITS;
                tmp[i + 3] += x >> 22;
                // At position 86, which is the starting bit position for word 3, we
                // have a factor of 0xffffc00 = 2**28 - 2**10
                if tmp[i + 3] < 0x10000000 {
                    set4 = 1;
                    tmp[i + 3] += 0x10000000 & x_mask;
                    tmp[i + 3] -= (x << 10) & LimbPattern::WIDTH28BITS;
                } else {
                    tmp[i + 3] -= (x << 10) & LimbPattern::WIDTH28BITS;
                }
                if tmp[i + 4] < 0x20000000 {
                    tmp[i + 4] += 0x20000000 & x_mask;
                    tmp[i + 4] -= set4;
                    tmp[i + 4] -= x >> 18;
                    if tmp[i + 5] < 0x10000000 {
                        tmp[i + 5] += 0x10000000 & x_mask;
                        tmp[i + 5] -= 1;
                        if tmp[i + 6] < 0x20000000 {
                            set7 = 1;
                            tmp[i + 6] += 0x20000000 & x_mask;
                            tmp[i + 6] -= 1;
                        } else {
                            tmp[i + 6] -= 1;
                        }
                    } else {
                        tmp[i + 5] -= 1;
                    }
                } else {
                    tmp[i + 4] -= set4;
                    tmp[i + 4] -= x >> 18;
                }
                // At position 200, which is the starting bit position for word 7, we
                // have a factor of 0xeffffff = 2**28 - 2**24 - 1
                if tmp[i + 7] < 0x10000000 {
                    tmp[i + 7] += 0x10000000 & x_mask;
                    tmp[i + 7] -= set7;
                    tmp[i + 7] -= (x << 24) & LimbPattern::WIDTH28BITS;
                    tmp[i + 8] += (x << 28) & LimbPattern::WIDTH29BITS;
                    if tmp[i + 8] < 0x20000000 {
                        tmp[i + 8] += 0x20000000 & x_mask;
                        tmp[i + 8] -= 1;
                        tmp[i + 8] -= x >> 4;
                        tmp[i + 9] += ((x >> 1) - 1) & x_mask;
                    } else {
                        tmp[i + 8] -= 1;
                        tmp[i + 8] -= x >> 4;
                        tmp[i + 9] += (x >> 1) & x_mask;
                    }
                } else {
                    tmp[i + 7] -= set7;
                    tmp[i + 7] -= (x << 24) & LimbPattern::WIDTH28BITS;
                    tmp[i + 8] += (x << 28) & LimbPattern::WIDTH29BITS;
                    if tmp[i + 8] < 0x20000000 {
                        tmp[i + 8] += 0x20000000 & x_mask;
                        tmp[i + 8] -= x >> 4;
                        tmp[i + 9] += ((x >> 1) - 1) & x_mask;
                    } else {
                        tmp[i + 8] -= x >> 4;
                        tmp[i + 9] += (x >> 1) & x_mask;
                    }
                }
            }

            if (i + 1) == 9 {
                break;
            }
            tmp[i + 2] += tmp[i + 1] >> 28;
            x = tmp[i + 1] & LimbPattern::WIDTH28BITS;
            tmp[i + 1] = 0;

            if x > 0 {
                let mut set5 = 0;
                let mut set8 = 0;
                let mut set9 = 0;
                x_mask = Self::mask(x);
                tmp[i + 3] += (x << 7) & LimbPattern::WIDTH28BITS;
                tmp[i + 4] += x >> 21;
                // At position 85, which is the starting bit position for word 3, we
                // have a factor of 0x1ffff800 = 2**29 - 2**11
                if tmp[i + 4] < 0x20000000 {
                    set5 = 1;
                    tmp[i + 4] += 0x20000000 & x_mask;
                    tmp[i + 4] -= (x << 11) & LimbPattern::WIDTH29BITS;
                } else {
                    tmp[i + 4] -= (x << 11) & LimbPattern::WIDTH29BITS;
                }
                if tmp[i + 5] < 0x10000000 {
                    tmp[i + 5] += 0x10000000 & x_mask;
                    tmp[i + 5] -= set5;
                    tmp[i + 5] -= x >> 18;
                    if tmp[i + 6] < 0x20000000 {
                        tmp[i + 6] += 0x20000000 & x_mask;
                        tmp[i + 6] -= 1;
                        if tmp[i + 7] < 0x10000000 {
                            set8 = 1;
                            tmp[i + 7] += 0x10000000 & x_mask;
                            tmp[i + 7] -= 1;
                        } else {
                            tmp[i + 7] -= 1;
                        }
                    } else {
                        tmp[i + 6] -= 1;
                    }
                } else {
                    tmp[i + 5] -= set5;
                    tmp[i + 5] -= x >> 18;
                }

                if tmp[i + 8] < 0x20000000 {
                    set9 = 1;
                    tmp[i + 8] += 0x20000000 & x_mask;
                    tmp[i + 8] -= set8;
                    tmp[i + 8] -= (x << 25) & LimbPattern::WIDTH29BITS;
                } else {
                    tmp[i + 8] -= set8;
                    tmp[i + 8] -= (x << 25) & LimbPattern::WIDTH29BITS;
                }
                if tmp[i + 9] < 0x10000000 {
                    tmp[i + 9] += 0x10000000 & x_mask;
                    tmp[i + 9] -= set9;
                    tmp[i + 9] -= x >> 4;
                    tmp[i + 10] += (x - 1) & x_mask;
                } else {
                    tmp[i + 9] -= set9;
                    tmp[i + 9] -= x >> 4;
                    tmp[i + 10] += x & x_mask;
                }
            }

            i += 2;
        }


        carry = 0;
        i = 0;
        while i < 8 {
            a.data[i] = tmp[i + 9];
            a.data[i] += carry;
            a.data[i] += (tmp[i + 10] << 28) & LimbPattern::WIDTH29BITS;
            carry = a.data[i] >> 29;
            a.data[i] &= LimbPattern::WIDTH29BITS;

            i += 1;
            a.data[i] = tmp[i + 9] >> 1;
            a.data[i] += carry;
            carry = a.data[i] >> 28;
            a.data[i] &= LimbPattern::WIDTH28BITS;

            i += 1;
        }

        a.data[8] = tmp[17];
        a.data[8] += carry;
        carry = a.data[8] >> 29;
        a.data[8] &= LimbPattern::WIDTH29BITS;

        Self::reduce_carry(a, carry as usize)
    }
}


impl Payload {
    fn init() -> Self {
        Payload { data: [0u32; 9] }
    }

    /// payload3 = payload1 + payload2
    ///
    /// payload1 = \[x0, x1, x2, x3, x4, x5, x6, x7, x8]
    /// payload2 = \[y0, y1, y2, y3, y4, y5, y6, y7, y8]
    ///
    /// |capacity| 32bits | 32bits | 32bits | 32bits | 32bits | 32bits | 32bits | 32bits | 32bits |
    /// |        |   x8   |   x7   |   x6   |   x5   |   x4   |   x3   |   x2   |   x1   |   x0   |
    /// |        |   y8   |   y7   |   y6   |   y5   |   y4   |   y3   |   y2   |   y1   |   y0   |
    ///
    /// |capacity| 32bits | 32bits | 32bits | 32bits | 32bits | 32bits | 32bits | 32bits | 32bits |
    /// | length | 29bits | 28bits | 29bits | 28bits | 29bits | 28bits | 29bits | 28bits | 29bits |
    /// |   257  |   228  |   200  |   171  |   143  |   114  |   86   |   57   |   29   |    0   |
    /// | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ |
    /// | carry  |   r8   |   r7   |   r6   |   r5   |   r4   |   r3   |   r2   |   r1   |   r0   |
    ///
    /// On entry, payload1\[i] + payload2\[i] must not overflow a 32-bit word.
    /// On exit: payload3\[0,2,...] < 2^30, payload3\[1,3,...] < 2^29
    fn add(&self, other: Payload) -> Payload {
        let mut result = Payload::init();
        let mut carry: u32 = 0;
        let mut i = 0;
        loop {
            let x = self.data[i].wrapping_add(other.data[i]).wrapping_add(carry);
            carry = x.shr(29);
            result.data[i] = x & (LimbPattern::WIDTH29BITS as u32);
            i += 1;
            if i == 9 {
                break;
            }
            let x = self.data[i].wrapping_add(other.data[i]).wrapping_add(carry);
            carry = x.shr(28);
            result.data[i] = x & (LimbPattern::WIDTH28BITS as u32);
            i += 1;
        }
        PayloadHelper::reduce_carry(&mut result, carry as usize);
        result
    }

    /// payload3 = payload1 - payload2
    ///
    /// On entry: payload1\[0,2,...] < 2^30, payload1\[1,3,...] < 2^29 and
    ///           payload2\[0,2,...] < 2^30, payload2\[1,3,...] < 2^29.
    /// On exit:  payload3\[0,2,...] < 2^30, payload3\[1,3,...] < 2^29.
    fn subtract(&self, other: Payload) -> Payload {
        let mut result = Payload::init();
        let mut carry: u32 = 0;
        let mut i = 0;
        loop {
            let x = self.data[i].wrapping_sub(other.data[i]).wrapping_add(P256ZERO31[i]).wrapping_add(carry);
            carry = x.shr(29);
            result.data[i] = x & (LimbPattern::WIDTH29BITS as u32);
            i += 1;
            if i == 9 {
                break;
            }
            let x = self.data[i].wrapping_sub(other.data[i]).wrapping_add(P256ZERO31[i]).wrapping_add(carry);
            carry = x.shr(28);
            result.data[i] = x & (LimbPattern::WIDTH28BITS as u32);
            i += 1;
        }
        PayloadHelper::reduce_carry(&mut result, carry as usize);
        result
    }

    /// multiply sets payload3 = payload1 * payload2.
    ///
    /// On entry: payload1\[0,2,...] < 2^30, payload1\[1,3,...] < 2^29 and
    ///           payload2\[0,2,...] < 2^30, payload2\[1,3,...] < 2^29.
    /// On exit:  payload3\[0,2,...] < 2^30, payload3\[1,3,...] < 2^29.
    fn multiply(&self, other: Payload) -> Payload {
        let mut result = Payload::init();
        let mut tmp: [u64; 17] = [0; 17];
        tmp[0] = (self.data[0] as u64) * (other.data[0] as u64);
        tmp[1] = (self.data[0] as u64) * ((other.data[1] as u64) << 0) +
            (self.data[1] as u64) * ((other.data[0] as u64) << 0);
        tmp[2] = (self.data[0] as u64) * ((other.data[2] as u64) << 0) +
            (self.data[1] as u64) * ((other.data[1] as u64) << 1) +
            (self.data[2] as u64) * ((other.data[0] as u64) << 0);
        tmp[3] = (self.data[0] as u64) * ((other.data[3] as u64) << 0) +
            (self.data[1] as u64) * ((other.data[2] as u64) << 0) +
            (self.data[2] as u64) * ((other.data[1] as u64) << 0) +
            (self.data[3] as u64) * ((other.data[0] as u64) << 0);
        tmp[4] = (self.data[0] as u64) * ((other.data[4] as u64) << 0) +
            (self.data[1] as u64) * ((other.data[3] as u64) << 1) +
            (self.data[2] as u64) * ((other.data[2] as u64) << 0) +
            (self.data[3] as u64) * ((other.data[1] as u64) << 1) +
            (self.data[4] as u64) * ((other.data[0] as u64) << 0);
        tmp[5] = (self.data[0] as u64) * ((other.data[5] as u64) << 0) +
            (self.data[1] as u64) * ((other.data[4] as u64) << 0) +
            (self.data[2] as u64) * ((other.data[3] as u64) << 0) +
            (self.data[3] as u64) * ((other.data[2] as u64) << 0) +
            (self.data[4] as u64) * ((other.data[1] as u64) << 0) +
            (self.data[5] as u64) * ((other.data[0] as u64) << 0);
        tmp[6] = (self.data[0] as u64) * ((other.data[6] as u64) << 0) +
            (self.data[1] as u64) * ((other.data[5] as u64) << 1) +
            (self.data[2] as u64) * ((other.data[4] as u64) << 0) +
            (self.data[3] as u64) * ((other.data[3] as u64) << 1) +
            (self.data[4] as u64) * ((other.data[2] as u64) << 0) +
            (self.data[5] as u64) * ((other.data[1] as u64) << 1) +
            (self.data[6] as u64) * ((other.data[0] as u64) << 0);
        tmp[7] = (self.data[0] as u64) * ((other.data[7] as u64) << 0) +
            (self.data[1] as u64) * ((other.data[6] as u64) << 0) +
            (self.data[2] as u64) * ((other.data[5] as u64) << 0) +
            (self.data[3] as u64) * ((other.data[4] as u64) << 0) +
            (self.data[4] as u64) * ((other.data[3] as u64) << 0) +
            (self.data[5] as u64) * ((other.data[2] as u64) << 0) +
            (self.data[6] as u64) * ((other.data[1] as u64) << 0) +
            (self.data[7] as u64) * ((other.data[0] as u64) << 0);
        tmp[8] = (self.data[0] as u64) * ((other.data[8] as u64) << 0) +
            (self.data[1] as u64) * ((other.data[7] as u64) << 1) +
            (self.data[2] as u64) * ((other.data[6] as u64) << 0) +
            (self.data[3] as u64) * ((other.data[5] as u64) << 1) +
            (self.data[4] as u64) * ((other.data[4] as u64) << 0) +
            (self.data[5] as u64) * ((other.data[3] as u64) << 1) +
            (self.data[6] as u64) * ((other.data[2] as u64) << 0) +
            (self.data[7] as u64) * ((other.data[1] as u64) << 1) +
            (self.data[8] as u64) * ((other.data[0] as u64) << 0);
        tmp[9] = (self.data[1] as u64) * ((other.data[8] as u64) << 0) +
            (self.data[2] as u64) * ((other.data[7] as u64) << 0) +
            (self.data[3] as u64) * ((other.data[6] as u64) << 0) +
            (self.data[4] as u64) * ((other.data[5] as u64) << 0) +
            (self.data[5] as u64) * ((other.data[4] as u64) << 0) +
            (self.data[6] as u64) * ((other.data[3] as u64) << 0) +
            (self.data[7] as u64) * ((other.data[2] as u64) << 0) +
            (self.data[8] as u64) * ((other.data[1] as u64) << 0);
        tmp[10] = (self.data[2] as u64) * ((other.data[8] as u64) << 0) +
            (self.data[3] as u64) * ((other.data[7] as u64) << 1) +
            (self.data[4] as u64) * ((other.data[6] as u64) << 0) +
            (self.data[5] as u64) * ((other.data[5] as u64) << 1) +
            (self.data[6] as u64) * ((other.data[4] as u64) << 0) +
            (self.data[7] as u64) * ((other.data[3] as u64) << 1) +
            (self.data[8] as u64) * ((other.data[2] as u64) << 0);
        tmp[11] = (self.data[3] as u64) * ((other.data[8] as u64) << 0) +
            (self.data[4] as u64) * ((other.data[7] as u64) << 0) +
            (self.data[5] as u64) * ((other.data[6] as u64) << 0) +
            (self.data[6] as u64) * ((other.data[5] as u64) << 0) +
            (self.data[7] as u64) * ((other.data[4] as u64) << 0) +
            (self.data[8] as u64) * ((other.data[3] as u64) << 0);
        tmp[12] = (self.data[4] as u64) * ((other.data[8] as u64) << 0) +
            (self.data[5] as u64) * ((other.data[7] as u64) << 1) +
            (self.data[6] as u64) * ((other.data[6] as u64) << 0) +
            (self.data[7] as u64) * ((other.data[5] as u64) << 1) +
            (self.data[8] as u64) * ((other.data[4] as u64) << 0);
        tmp[13] = (self.data[5] as u64) * ((other.data[8] as u64) << 0) +
            (self.data[6] as u64) * ((other.data[7] as u64) << 0) +
            (self.data[7] as u64) * ((other.data[6] as u64) << 0) +
            (self.data[8] as u64) * ((other.data[5] as u64) << 0);
        tmp[14] = (self.data[6] as u64) * ((other.data[8] as u64) << 0) +
            (self.data[7] as u64) * ((other.data[7] as u64) << 1) +
            (self.data[8] as u64) * ((other.data[6] as u64) << 0);
        tmp[15] = (self.data[7] as u64) * ((other.data[8] as u64) << 0) +
            (self.data[8] as u64) * ((b[7] as u64) << 0);
        tmp[16] = (self.data[8] as u64) * ((other.data[8] as u64) << 0);

        PayloadHelper::reduce_degree(&mut result, &mut tmp);
        result
    }

    fn square(&self) -> Payload {
        let mut result = Payload::init();
        let mut tmp: [u64; 17] = [0; 17];
        tmp[0] = (self.data[0] as u64) * (self.data[0] as u64);
        tmp[1] = (self.data[0] as u64) * ((self.data[1] as u64) << 1);
        tmp[2] = (self.data[0] as u64) * ((self.data[2] as u64) << 1) +
            (self.data[1] as u64) * ((self.data[1] as u64) << 1);
        tmp[3] = (self.data[0] as u64) * ((self.data[3] as u64) << 1) +
            (self.data[1] as u64) * ((self.data[2] as u64) << 1);
        tmp[4] = (self.data[0] as u64) * ((self.data[4] as u64) << 1) +
            (self.data[1] as u64) * ((self.data[3] as u64) << 2) +
            (self.data[2] as u64) * (self.data[2] as u64);
        tmp[5] = (self.data[0] as u64) * ((self.data[5] as u64) << 1) +
            (self.data[1] as u64) * ((self.data[4] as u64) << 1) +
            (self.data[2] as u64) * ((self.data[3] as u64) << 1);
        tmp[6] = (self.data[0] as u64) * ((self.data[6] as u64) << 1) +
            (self.data[1] as u64) * ((self.data[5] as u64) << 2) +
            (self.data[2] as u64) * ((self.data[4] as u64) << 1) +
            (self.data[3] as u64) * ((self.data[3] as u64) << 1);
        tmp[7] = (self.data[0] as u64) * ((self.data[7] as u64) << 1) +
            (self.data[1] as u64) * ((self.data[6] as u64) << 1) +
            (self.data[2] as u64) * ((self.data[5] as u64) << 1) +
            (self.data[3] as u64) * ((self.data[4] as u64) << 1);
        tmp[8] = (self.data[0] as u64) * ((self.data[8] as u64) << 1) +
            (self.data[1] as u64) * ((self.data[7] as u64) << 2) +
            (self.data[2] as u64) * ((self.data[6] as u64) << 1) +
            (self.data[3] as u64) * ((self.data[5] as u64) << 2) +
            (self.data[4] as u64) * (self.data[4] as u64);
        tmp[9] = (self.data[1] as u64) * ((self.data[8] as u64) << 1) +
            (self.data[2] as u64) * ((self.data[7] as u64) << 1) +
            (self.data[3] as u64) * ((self.data[6] as u64) << 1) +
            (self.data[4] as u64) * ((self.data[5] as u64) << 1);
        tmp[10] = (self.data[2] as u64) * ((self.data[8] as u64) << 1) +
            (self.data[3] as u64) * ((self.data[7] as u64) << 2) +
            (self.data[4] as u64) * ((self.data[6] as u64) << 1) +
            (self.data[5] as u64) * ((self.data[5] as u64) << 1);
        tmp[11] = (self.data[3] as u64) * ((self.data[8] as u64) << 1) +
            (self.data[4] as u64) * ((self.data[7] as u64) << 1) +
            (self.data[5] as u64) * ((self.data[6] as u64) << 1);
        tmp[12] = (self.data[4] as u64) * ((self.data[8] as u64) << 1) +
            (self.data[5] as u64) * ((self.data[7] as u64) << 2) +
            (self.data[6] as u64) * (self.data[6] as u64);
        tmp[13] = (self.data[5] as u64) * ((self.data[8] as u64) << 1) +
            (self.data[6] as u64) * ((self.data[7] as u64) << 1);
        tmp[14] = (self.data[6] as u64) * ((self.data[8] as u64) << 1) +
            (self.data[7] as u64) * ((self.data[7] as u64) << 1);
        tmp[15] = (self.data[7] as u64) * ((self.data[8] as u64) << 1);
        tmp[16] = (self.data[8] as u64) * (self.data[8] as u64);

        PayloadHelper::reduce_degree(&mut result, &mut tmp);
        result
    }
}

#[cfg(test)]
mod tests {
    use std::ops::Sub;
    use num_bigint::Sign;
    use num_traits::Num;

    use super::*;

    #[test]
    fn ri() {
        let r_hex = "20000000000000000000000000000000000000000000000000000000000000000";
        let ri_hex = "7ffffffd80000002fffffffe000000017ffffffe800000037ffffffc80000002";

        let p = BigInt::from_bytes_be(Sign::Plus, &EC_P);
        let r = BigInt::from(2u64).pow(257);
        assert_eq!(r.to_str_radix(16), r_hex);

        let ext = r.extended_gcd(&p);
        let ri = ext.x.mod_floor(&p);
        assert_eq!(ri.to_str_radix(16), ri_hex);
    }

    #[test]
    fn payload() {
        let n = "115792089210356248756420345214020892766250353991924191454421193933289684991996";
        let n = BigInt::from_str_radix(n, 10).unwrap();
        let payload = PayloadHelper::transform(&n);
        assert_eq!(payload.data, [
            536870905, 268435455, 895,
            268428288, 536870911, 268435455,
            536870911, 150994943, 268435455
        ]);
        let m = PayloadHelper::restore(&payload);
        assert_eq!(m, n);
    }

    #[test]
    fn carry_table() {
        let mut table: [u32; 8 * 9] = [0; 72];
        for i in 0..8 {
            let value = BigInt::from(i as i64);
            let payload = PayloadHelper::transform(&value);
            for (j, e) in payload.data.iter().enumerate() {
                table[i * 9 + j] = *e;
                print!("0x{:>08X}, ", *e);
            }
            println!();
        }
        assert_eq!(table, P256CARRY)
    }

    #[test]
    fn factor_table() {
        let mut table: [[u32; 9]; 9] = [[0; 9]; 9];
        for i in 0..9 {
            let value = BigInt::from(i as i64);
            let payload = PayloadHelper::transform(&value);

            let mut temp: [u32; 9] = [0; 9];
            for (j, e) in payload.data.iter().enumerate() {
                temp[j] = *e;
                print!("0x{:>08X}, ", *e);
            }
            table[i] = temp;
            println!();
        }
        assert_eq!(table, P256FACTOR)
    }

    #[test]
    fn zero31() {
        let limbs_to_big = |p: [u32; 9]| {
            let mut n = BigInt::from_u32(p[8]).unwrap();
            let mut i: isize = 7;
            while i >= 0 {
                if (i & 1) == 0 {
                    n = n.shl(29);
                } else {
                    n = n.shl(28);
                }
                n = n.add(BigInt::from_u32(p[i as usize]).unwrap());
                i -= 1;
            }
            n
        };
        let result = limbs_to_big(
            [1 << 31, 1 << 30, 1 << 31, 1 << 30, 1 << 31, 1 << 30, 1 << 31, 1 << 30, 1 << 31]
        );
        let p256 = P256Elliptic::init();
        let mut temp = result.mod_floor(&p256.ec.p.to_bigint().unwrap());
        temp = result.sub(temp);

        let mut out: [u32; 9] = [0; 9];
        let mut i: usize = 0;
        while i < 9 {
            let bits = temp.to_u64_digits().1;
            if !bits.is_empty() {
                out[i] = (bits[0] as u32) & 0x7fffffff;
                if out[i] < 0x70000000 {
                    out[i] += 0x80000000;
                }
            } else {
                out[i] = 0x80000000
            }
            temp = temp.sub(BigInt::from_u64(out[i] as u64).unwrap());
            temp = temp.shr(29);
            i += 1;
            if i == 9 {
                break;
            }
            let bits = temp.to_u64_digits().1;
            if !bits.is_empty() {
                out[i] = (bits[0] as u32) & 0x3fffffff;
                if out[i] < 0x30000000 {
                    out[i] += 0x40000000;
                }
            } else {
                out[i] = 0x40000000;
            }
            temp = temp.sub(BigInt::from_u64(out[i] as u64).unwrap());
            temp = temp.shr(28);
            i += 1;
        }

        println!("{:08X?}", out);
        assert_eq!(out, P256ZERO31);
    }
}