# SPDX-License-Identifier: MIT
#
# Generate response KAT digest from SageMath SNOVA
#
# Copyright (c) 2025 SNOVA TEAM

import math
import os
import subprocess
import sys

from hashlib import shake_256


variants = [
    [28, 5, 19, 4],
    [24, 5, 16, 4],
    [48, 17, 16, 2],
    [48, 16, 19, 2],
    [28, 4, 16, 4, 5],
    [28, 4, 19, 4, 5],

    [40, 7, 19, 4],
    [37, 8, 16, 4],
    [72, 25, 16, 2],
    [72, 24, 19, 2],
    [38, 5, 16, 4, 5],
    [38, 5, 19, 4, 5],

    [50, 9, 19, 4],
    [60, 10, 16, 4],
    [97, 33, 16, 2],
    [96, 32, 19, 2],
    [52, 6, 16, 4, 6],
    [52, 6, 19, 4, 6],
]


for var in variants:
    v = var[0]
    o = var[1]
    q = var[2]
    l = var[3]

    if len(var) > 4:
        r = var[4]
    else:
        r = l

    if len(var) > 5:
        m1 = var[5]
    else:
        m1 = math.ceil(o * r / l)

    if len(var) > 6:
        n_alpha = var[6]
    else:
        n_alpha = l * r + r

    for aes in [True, False]:
        name = f"SNOVA_{v}_{o}_{q}_{l if l == r else f'{l}x{r}'}{'_AES' if aes else ''}"

        with open('snova.sage') as infile:
            data = infile.read()

        # Create sage file from parameters
        param_dict = {
            "print('# SNOVA', v, o, q, l, 'AES' if aes else 'SHAKE', r, m1, n_alpha)":
            f"print('# {name}')",
            "v = 28": f"v = {v}",
            "o = 5": f"o = {o}",
            "q = 19": f"q = {q}",
            "l = 4": f"l = {l}",
            "aes = False": f"aes = {aes}",
            "r = l": f"r = {r}",
            "m1 = math.ceil(o * r / l)": f"m1 = {m1}",
            "n_alpha = r * r + r": f"n_alpha = {n_alpha}",
            "range(1)": "range(100)",
        }
        for key in param_dict.keys():
            data = data.replace(key, param_dict[key])

        with open(f'_{name}.sage', 'w') as outfile:
            outfile.write(data)

        # Run and print digest

        result = subprocess.run(['sage', f'_{name}.sage'], capture_output=True, text=True)
        digest = shake_256(result.stdout.encode()).digest(24).hex()
        print(f'{name}.rsp  {digest}')
        sys.stdout.flush()

        os.remove(f'_{name}.sage')
