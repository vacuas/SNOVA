import shutil
import os

params = [
    [24, 5, 23, 4, True],
    [24, 5, 23, 4, False],
    [24, 5, 16, 4, True],
    [24, 5, 16, 4, False],
    [36, 7, 19, 2, True, 5, 16, 20],
    [36, 7, 19, 2, False, 5, 16, 20],
    [43, 17, 16, 2, True],
    [43, 17, 16, 2, False],

    [37, 8, 19, 4, True],
    [37, 8, 19, 4, False],
    [37, 8, 16, 4, True],
    [37, 8, 16, 4, False],
    [52, 8, 19, 2, True, 6, 24, 24],
    [52, 8, 19, 2, False, 6, 24, 24],
    [69, 25, 16, 2, True],
    [69, 25, 16, 2, False],

    [60, 10, 23, 4, True],
    [60, 10, 23, 4, False],
    [60, 10, 16, 4, True],
    [60, 10, 16, 4, False],
    [70, 11, 19, 2, True, 6, 32, 30],
    [70, 11, 19, 2, False, 6, 32, 30],
    [99, 33, 16, 2, True],
    [99, 33, 16, 2, False],
]


mf = open('Makefile.root', 'w')
print('''MAKEFLAGS += --no-print-directory

all: kat

digest:
	@sh digest.sh */PQCsign*''', file=mf)

print('\nkat:', file=mf)
for param in params:
    v = param[0]
    o = param[1]
    q = param[2]
    l = param[3]
    aes = param[4]
    print('''\tmake -C SNOVA_{}_{}_{}_{}{} kat'''.format(
        v, o, q, l, '_AES' if aes else ''), file=mf)

print('\nspeed:', file=mf)
for param in params:
    v = param[0]
    o = param[1]
    q = param[2]
    l = param[3]
    aes = param[4]
    print('''\t@make -C SNOVA_{}_{}_{}_{}{} speed'''.format(
        v, o, q, l, '_AES' if aes else ''), file=mf)

print('\nclean:', file=mf)
for param in params:
    v = param[0]
    o = param[1]
    q = param[2]
    l = param[3]
    aes = param[4]
    print('''\tmake -C SNOVA_{}_{}_{}_{}{} clean'''.format(
        v, o, q, l, '_AES' if aes else ''), file=mf)

mf.close()

ref_sources = [
    '.gitignore',
    'LICENSE',
    'PQCgenKAT_sign.c',
    'aes.c',
    'api.h',
    'keccak_opt64.h',
    'rng.c',
    'rng.h',
    'sign.c',
    'snova.h',
    'snova_ref.c',
    'speed.c',
    'symmetric.h',
    'symmetric_ref.c',
]

source_dir = '../src/'

for target in ['ref', 'opt', 'avx2']:
    shutil.rmtree(target, ignore_errors=True)
    os.makedirs(target)
    shutil.copy('Makefile.root', target + '/Makefile')
    shutil.copy('digest.sh', target)

    for param in params:
        v = param[0]
        o = param[1]
        q = param[2]
        l = param[3]
        aes = param[4]

        dirname = target + '/SNOVA_{}_{}_{}_{}{}/'.format(
            v, o, q, l, '_AES' if aes else '')

        if target == 'ref':
            os.makedirs(dirname)
            shutil.copy('README.ref', dirname + 'README.md')
            shutil.copy('Makefile.ref', dirname + 'Makefile')
            for file in ref_sources:
                shutil.copy(source_dir + file, dirname)
        else:
            shutil.copytree(source_dir, dirname)

        if target == 'opt':
            os.remove(dirname + 'snova_avx2.c')
            os.remove(dirname + 'snova_avx2_q.c')
            os.remove(dirname + 'snova_avx2_16.c')
            os.remove(dirname + 'snova_avx2_rect.c')
            shutil.copy('README.opt', dirname + 'README.md')
            shutil.copy('Makefile.opt', dirname + 'Makefile')

        sp = open(dirname + 'snova_params.h', 'w')
        print('#define SNOVA_v', v, file=sp)
        print('#define SNOVA_o', o, file=sp)
        print('#define SNOVA_q', q, file=sp)
        print('#define SNOVA_l', l, file=sp)
        if aes:
            print('#define AESCTR', file=sp)
        if len(param) > 5:
            print('#define SNOVA_r', param[5], file=sp)
        if len(param) > 6:
            print('#define SNOVA_m1', param[6], file=sp)
        if len(param) > 7:
            print('#define SNOVA_alpha', param[7], file=sp)
        sp.close()

os.remove('Makefile.root')
