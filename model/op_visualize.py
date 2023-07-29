# some codes are from official tutorials https://captum.ai/tutorials/Segmentation_Interpret

import os
import cv2 as cv
import numpy as np
import torch
import shutil
import utils.image_transforms as image_transforms
from torch.utils.data import DataLoader
import op_dataset as dataset
from utils.loss_function import *
from model_unet import Unet
from utils import image_helpers
from PIL import Image
import os
from utils import image_helpers
import os
import numpy as np
from captum.attr import visualization as viz
from captum.attr import LayerGradCam, FeatureAblation, LayerActivation, LayerAttribution, NoiseTunnel
from model_fcn import VGGNet, FCNs
from model_deeplabv3 import DeepLabV3
from model_liteunet import LiteUnet
from captum.attr import IntegratedGradients
from matplotlib.colors import LinearSegmentedColormap
from model_miniunet import MiniUnet

os.environ["CUDA_VISIBLE_DEVICES"] = "3"


def visualize_via_captum(data_path, model_path, crop_size, device="cpu"):

    def agg_segmentation_wrapper(inp):
        model_out = model(inp)
        # Creates binary matrix with 1 for original argmax class for each pixel
        # and 0 otherwise. Note that this may change when the input is ablated
        # so we use the original argmax predicted above, out_max.
        selected_inds = torch.zeros_like(model_out[0:1]).scatter_(1, out_max, 1)
        return (model_out * selected_inds).sum(dim=(2, 3))

    # # STEP: load data
    image_norm_transform = None
    center_crop_transform = image_transforms.CenterCrop(crop_size)
    img2tensor_transform = image_transforms.NpyToTensor()
    mask2tensor_transform = image_transforms.MaskToTensor()

    visual_set = dataset.Dataset(data_path, 'visual', 1, image_norm_transform=image_norm_transform, img2tensor_transform=img2tensor_transform, center_crop_transform=center_crop_transform, mask2tensor_transform=mask2tensor_transform)
    visual_loader = DataLoader(visual_set, batch_size=1, shuffle=False)

    palette = dataset.palette
    num_classes = dataset.num_classes

    # # STEP: load model
    # model = MiniUnet(img_ch=3, num_classes=num_classes, depth=2)

    model = LiteUnet(img_ch=3, num_classes=num_classes, depth=2)
    #
    # model = DeepLabV3(num_classes=dataset.num_classes)

    # vgg_model = VGGNet(requires_grad=True, remove_fc=True)
    # model = FCNs(pretrained_net=vgg_model, n_class=dataset.num_classes)

    model.load_state_dict(torch.load(model_path, map_location=torch.device('cpu')))
    model.eval()

    for (input, mask), file_name in visual_loader:
        file_name = file_name[0].split('/')[-1]
        print(file_name)
        if device == "gpu":
            X = input.cuda()
        else:
            X = input

        model_out = model(X)
        # model_out = torch.sigmoid(model_out)

        mask_matrix = image_helpers.onehot_to_mask(np.array(mask).squeeze().transpose([1, 2, 0]), palette)[:, :, 1]
        mask_matrix[mask_matrix == 255] = 0
        mask_matrix[mask_matrix == 146] = 1
        mask_matrix[mask_matrix == 246] = 2
        mask_matrix[mask_matrix == 229] = 3
        mask_matrix[mask_matrix == 106] = 4
        mask_matrix[mask_matrix == 150] = 5
        mask_matrix[mask_matrix == 50] = 6

        out_max = torch.argmax(model_out, dim=1, keepdim=True)

        sem_classes = ['BG', 'REF', 'INS', 'DEL', 'INV', 'DUP', 'invDUP']

        # # STEP: visualize layer grad cam
        target_layver = model.conv_1x1

        target_class = 3

        """ old color
        # default_cmap = LinearSegmentedColormap.from_list('custom', [(0, '#fff3f2'), (0.25, '#feebd6'), (0.85, '#ff9f9f'), (1, '#9e0021')], N=256)
        # 
        # default_cmap = LinearSegmentedColormap.from_list('custom', [(0, '#fff7ec'), (0.25, '#fdd09a'), (0.85, '#e04630'), (1, '#800000')], N=256)
        # default_cmap = LinearSegmentedColormap.from_list('custom', [(0, '#fff7ec'), (0.25, '#fff7ec'), (0.75, '#e04630'), (1, '#800000')], N=256)
        # 
        # default_cmap = LinearSegmentedColormap.from_list('custom', [(0, '#fff5f1'), (0.25, '#ede9d0'), (0.5, '#cfa195'), (0.75, '#ff918d'), (1, '#ff8066')], N=256)
        # default_cmap = LinearSegmentedColormap.from_list('custom', [(0, '#ffffff'), (0.25, '#fff5f1'), (0.35, '#ede9d0'), (0.65, '#cfa195'), (0.75, '#ff918d'),
        # (1, '#ff8066')], N=256)
        """

        default_cmap = LinearSegmentedColormap.from_list('custom', [(0, '#fff7ec'), (0.4, '#fdd09a'), (0.8, '#e04630'), (1, '#800000')], N=256)

        lgc = LayerGradCam(agg_segmentation_wrapper, target_layver)
        gc_attr = lgc.attribute(input, target=target_class)

        lgc_img = viz.visualize_image_attr(gc_attr[0].cpu().permute(1, 2, 0).detach().numpy(), sign="all", cmap=default_cmap, show_colorbar=True)
        # lgc_img = viz.visualize_image_attr(gc_attr[0].cpu().permute(1, 2, 0).detach().numpy(), sign="all", show_colorbar=True)

        lgc_img[0].savefig("captum_lgc.png")


        # # STEP: do feature ablation
        fa = FeatureAblation(agg_segmentation_wrapper)
        fa_attr = fa.attribute(input, feature_mask=out_max, target=target_class)
        # fa_attr = fa.attribute(input, feature_mask=torch.tensor(mask_matrix), perturbations_per_eval=2, target=1)

        # fa_img = viz.visualize_image_attr(fa_attr[0].cpu().detach().permute(1, 2, 0).numpy(), sign="all", cmap=default_cmap, show_colorbar=True)
        fa_img = viz.visualize_image_attr(fa_attr[0].cpu().detach().permute(1, 2, 0).numpy(), sign="all", show_colorbar=True)
        fa_img[0].savefig("captum_fa.png")
        #
        #
        # # # STEP: visualize integrated gradients
        # integrated_gradients = IntegratedGradients(agg_segmentation_wrapper)
        # # ig_attr = integrated_gradients.attribute(input, target=1, baselines=0,  n_steps=200)
        # ig_attr = integrated_gradients.attribute(input, target=target_class, baselines=0)
        #
        # ig_img = viz.visualize_image_attr(ig_attr[0].cpu().permute(1, 2, 0).detach().numpy(), cmap=default_cmap, show_colorbar=True)
        # # ig_img = viz.visualize_image_attr(ig_attr[0].cpu().permute(1, 2, 0).detach().numpy(), sign="all",  show_colorbar=True)
        #
        # ig_img[0].savefig("captum_ig.png")


        # # STEP: visualize layer activation
        # la = LayerActivation(agg_segmentation_wrapper, target_layver)
        # activation = la.attribute(input)
        # la_img = viz.visualize_image_attr(activation[0][1:2, :, :].cpu().permute(1, 2, 0).detach().numpy(), sign="all", cmap=default_cmap,show_colorbar=True)
        # la_img[0].savefig("captum_la.png")


        # # STEP: smooth integrated gradients
        # noise_tunnel = NoiseTunnel(integrated_gradients)
        # nt_attr = noise_tunnel.attribute(input, nt_samples=10, nt_type='smoothgrad_sq', target=1)
        # nt_img = viz.visualize_image_attr(nt_attr[0].cpu().permute(1, 2, 0).detach().numpy(), sign="all",  show_colorbar=True)
        #
        # nt_img[0].savefig("captum_nt.png")

        # # STEP: visualize via occlusion
        # from captum.attr import Occlusion   # successful in node3, but cost a lot of memories and cpus
        # occlusion = Occlusion(model)
        # occ_attr = occlusion.attribute(input, strides=(3, 8, 8), target=1, sliding_window_shapes=(3, 15, 15), baselines=0)
        # occ_img = viz.visualize_image_attr(occ_attr[0].cpu().permute(1, 2, 0).detach().numpy(), sign="all",  show_colorbar=True)
        # occ_img[0].savefig("captum_occ.png")


if __name__ == '__main__':
    data_path = "data"

    # model_path = "/mnt/c/workspace/test/sim_train_ccs/model/mininet_batch1/model_miniunet_256_8_16_32.pth"

    # model_path = "/mnt/c/workspace/test/sim_simple_train_ccs/model/batch1/model_miniunet_256_8_16_32.pth"

    # model_path = "/data/home/songbo/workspace/svision-pro/sim_train_ccs/models/model_miniunet_256_8_16_32.pth"

    image_size = 1024

    model_path = "/mnt/c/workspace/test/sim_simple_train_ccs/models/model_liteunet_{}_8_16_32_32_32.pth".format(image_size)
    # model_path = "/data/home/songbo/workspace/svision-pro/sim_simple_train_ccs/models/model_liteunet_{}_8_16_32_32_32.pth".format(image_size)

    device = "cpu"

    visualize_via_captum(data_path, model_path, image_size, device)
    # visualize_via_gradcam(data_path, model_path, image_size, device)