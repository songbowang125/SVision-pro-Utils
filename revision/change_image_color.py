import os
import cv2 as cv
import numpy as np

if __name__ == '__main__':
    raw_img = "/mnt/d/workspace/svision-pro/sim_fig1/chr1_labels/sample.chr1_25557927_25561484_NA_4438_base0.png"

    out_img = raw_img.replace(".png", ".changed.png")


    raw_img_rgb = cv.imread(raw_img)

    # print(raw_img_rgb)
    print(np.shape(raw_img_rgb))

    raw_inv = np.where(raw_img_rgb==(187, 106, 127))
    for x, y in zip(raw_inv[0], raw_inv[1]):
        # raw_img_rgb[x, y, :] = (127, 127, 0)
        raw_img_rgb[x, y, :] = (89, 121, 204)

    raw_wt = np.where(raw_img_rgb==(197, 146, 66))
    for x, y in zip(raw_wt[0], raw_wt[1]):
        # raw_img_rgb[x, y, :] = (151, 84, 48)
        # raw_img_rgb[x, y, :] = (89, 89, 89)
        raw_img_rgb[x, y, :] = (127, 127, 127)

    cv.imwrite(out_img, raw_img_rgb)
    # for x, y, z in np.where(raw_img_rgb==(127, 106, 187)):
    #
    #     print(raw_img_rgb[x, y, :])