from .aes256_ctr_drbg import AES256_CTR_DRBG
from . import pyaes

def rng(seed):
    return AES256_CTR_DRBG(bytes(seed))

def aesctr(seed, num):
    res = bytearray()
    block_i = 0
    cipher = pyaes.AESModeOfOperationECB(bytes(seed))
    while len(res) < num:
        blockseed = bytearray()
        for j in range(16):
            blockseed.append((block_i >> (8 * (15 - j))) % 256)
        res += cipher.encrypt(bytes(blockseed))
        block_i += 1
    return bytes(res[:num])
