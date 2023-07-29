
from torchsummary import summary

import op_dataset as dataset
from model_unet import Unet
from model_liteunet import LiteUnet

from model_fcn import VGGNet, FCNs
from model_deeplabv3 import DeepLabV3
from model_miniunet import MiniUnet
from model_alexnet import AlexNet

def model_summary(num_classes =dataset.num_classes):
    # model = LiteUnet(img_ch=3, num_classes=num_classes, depth=2)
    model = AlexNet()

    # model = Unet(img_ch=3, num_classes=num_classes, depth=2)

    # model = DeepLabV3(num_classes=dataset.num_classes)

    # vgg_model = VGGNet(requires_grad=True, remove_fc=True)
    # model = FCNs(pretrained_net=vgg_model, n_class=dataset.num_classes)

    # summary(model, (3, 256, 256))
    summary(model, (3, 224, 224))

if __name__ == '__main__':


    model_summary()
