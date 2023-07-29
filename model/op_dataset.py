import os
from sklearn.model_selection import KFold
import os
import numpy as np
from PIL import Image
from torch.utils import data
from utils import image_helpers


class_colors = {
    "BG": [255, 255, 255],
    "REF": [66, 146, 197],
    "INS": [251, 246, 180],
    "DEL": [253, 229, 217],
    "INV": [127, 106, 187],
    "DUP": [235, 150, 0],
    "invDUP": [235, 50, 0]
}


#
palette = [class_colors["BG"], class_colors["REF"], class_colors["INS"], class_colors["DEL"], class_colors["INV"], class_colors["DUP"], class_colors["invDUP"]]  # one-hot
num_classes = len(palette)


def dataset_kfold(dataset_dir, save_path, k=5):
    """
    split dataset to perform kfold
    """

    data_list = os.listdir(dataset_dir)

    kf = KFold(k, shuffle=False)

    for i, (tr, val) in enumerate(kf.split(data_list), 1):
        print(len(tr), len(val))

        # # STEP: remove old files
        if os.path.exists(os.path.join(save_path, 'model{}.txt'.format(i))):
            os.remove(os.path.join(save_path, 'model{}.txt'.format(i)))
            os.remove(os.path.join(save_path, 'val{}.txt'.format(i)))

        # # STEP: output to files
        for item in tr:
            file_name = data_list[item]
            with open(os.path.join(save_path, 'model{}.txt'.format(i)), 'a') as f:
                f.write(file_name)
                f.write('\n')

        for item in val:
            file_name = data_list[item]
            with open(os.path.join(save_path, 'val{}.txt'.format(i)), 'a') as f:
                f.write(file_name)
                f.write('\n')


def make_dataset(root, mode, fold):
    assert mode in ['model', 'val', 'test', 'visual']

    items = []
    if mode == 'model':
        img_path = os.path.join(root, 'Images')
        mask_path = os.path.join(root, 'Labels')

        if 'Augdata' in root:  # 当使用增广后的训练集
            data_list = os.listdir(os.path.join(root, 'Labels'))
        else:
            data_list = [l.strip('\n') for l in open(os.path.join(root, 'model{}.txt'.format(fold))).readlines()]
        for it in data_list:
            item = (os.path.join(img_path, it), os.path.join(mask_path, it))
            items.append(item)

    elif mode == 'val':
        img_path = os.path.join(root, 'Images')
        mask_path = os.path.join(root, 'Labels')
        data_list = [l.strip('\n') for l in open(os.path.join(root, 'val{}.txt'.format(fold))).readlines()]
        for it in data_list:
            item = (os.path.join(img_path, it), os.path.join(mask_path, it))
            items.append(item)

    elif mode == "visual":
        img_path = os.path.join(root, 'Images')
        mask_path = os.path.join(root, 'Labels')
        data_list = [l.strip('\n') for l in open(os.path.join(root, 'visual{}.txt'.format(fold))).readlines()]
        for it in data_list:
            item = (os.path.join(img_path, it), os.path.join(mask_path, it))
            items.append(item)
    else:
        img_path = os.path.join(root, 'Images')
        mask_path = os.path.join(root, 'Labels')

        data_list = [l.strip('\n') for l in open(os.path.join(root, 'test.txt')).readlines()]
        for it in data_list:
            # item = (os.path.join(img_path, 'c0', it))
            item = (os.path.join(img_path, it), os.path.join(img_path, it))
            items.append(item)

    return items


class Dataset(data.Dataset):
    def __init__(self, root, mode, fold, image_norm_transform=None, center_crop_transform=None, img2tensor_transform=None, mask2tensor_transform=None):
        self.imgs = make_dataset(root, mode, fold)

        self.palette = palette
        self.mode = mode
        if len(self.imgs) == 0:
            raise RuntimeError('Found 0 images, please check the data set')
        self.mode = mode
        self.image_norm_transform = image_norm_transform
        self.center_crop_transform = center_crop_transform
        self.img2tensor_transform = img2tensor_transform
        self.mask2tensor_transform = mask2tensor_transform

    def __getitem__(self, index):

        img_path, mask_path = self.imgs[index]
        # print(img_path)
        file_name = mask_path.split('\\')[-1]

        # # load and precess img
        img = Image.open(img_path)
        if self.center_crop_transform is not None:
            img = self.center_crop_transform(img)

        img = np.array(img)
        img = img.transpose([2, 0, 1])

        if self.img2tensor_transform is not None:
            img = self.img2tensor_transform(img)

        if self.image_norm_transform is not None:
            img = self.image_norm_transform(img)

        # # load and precess mask
        if self.mode == "test":
            return img, file_name

        else:
            mask = Image.open(mask_path)
            if self.center_crop_transform is not None:
                mask = self.center_crop_transform(mask)

            mask = np.array(mask)

            mask = image_helpers.mask_to_onehot(mask, self.palette)

            # shape from (H, W, C) to (C, H, W)
            mask = mask.transpose([2, 0, 1])

            if self.mask2tensor_transform is not None:
                mask = self.mask2tensor_transform(mask)

            return (img, mask), file_name

    def __len__(self):
        return len(self.imgs)


def re_group_images(img_path, out_path):

    out_image_path = os.path.join(out_path, "Images")
    out_label_path = os.path.join(out_path, "Labels")

    if not os.path.exists(out_image_path):
        os.mkdir(out_image_path)

    if not os.path.exists(out_label_path):
        os.mkdir(out_label_path)

    for sub_interval_folder in os.listdir(img_path):

        origin_img_path = os.path.join(img_path, "{}/Images".format(sub_interval_folder ))

        for img in os.listdir(origin_img_path):

            if "segmentation" in img:
                os.system("mv {} {}".format(os.path.join(origin_img_path, img), os.path.join(out_label_path, img.replace(".segmentation", ""))))
            else:
                os.system("mv {} {}".format(os.path.join(origin_img_path, img), os.path.join(out_image_path, img)))


if __name__ == '__main__':


    data_path = "/data/home/songbo/workspace/svision-pro/sim_simple_ccs/sv_call"

    re_group_images("/data/home/songbo/workspace/svision-pro/sim_simple_ccs/sv_call/tmp_sim_221107_161737", data_path)

    # # split dataset to kfold
    dataset_kfold(os.path.join(data_path, 'Labels'), data_path)
