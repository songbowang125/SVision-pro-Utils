import time
import os
import torch
import random
from torch.utils.data import DataLoader
from torch.optim import lr_scheduler
from tqdm import tqdm
import sys
import datetime
import numpy as np
import op_dataset as dataset
import utils.image_transforms as image_transforms
from utils.loss_function import Dice_Loss, CE_Dice_Loss, multi_class_dice
from utils.early_stopping import EarlyStopping
from utils.LRScheduler import PolyLR
import os
from model_unet import Unet
import torch.nn as nn
import torch.optim as optim
from model_fcn import VGGNet, FCNs
from model_deeplabv3 import DeepLabV3
from model_miniunet import MiniUnet
from model_liteunet import LiteUnet
import torch.nn.functional as F

"""
    Setting hyper parameters
"""
crop_size = 1024  # # img crop size

output_dice = True

train_path = "/data/home/songbo/workspace/svision-pro/sim_simple_train_ccs/images_{}".format(crop_size)

val_path=train_path

device = "gpu"

model = "unet"

fold = 1

batch_size = 2  # batch size

n_epoch = 300  # max training epochs

os.environ["CUDA_VISIBLE_DEVICES"] = "3"

early_stop__eps = 1e-3  # early stopping threshhold
early_stop_patience = 15  # early stopping waits

initial_lr = 1e-4  # initial learning rate
threshold_lr = 1e-6  # learning rate for early stopping
weight_decay = 1e-5  # learning rate decay

model_name = '{}_fold{}_batch{}_{}'.format(model, fold, batch_size,  random.randint(1, 1e6))

log_file = open(os.path.join(train_path, "{}.txt".format(model_name)), "w")


"""
    End
"""


def train(train_loader, val_loader, net, criterion, optimizer, scheduler, early_stopping, num_epoches, iters):

    for epoch in range(1, num_epoches + 1):
        st = datetime.datetime.now()

        train_losses = []
        val_losses = []

        # model model
        net.train().cuda()

        epoch_classes_dice = [[] for i in range(dataset.num_classes)]
        epoch_mean_dice = 0

        for batch, ((input, mask), file_name) in enumerate(train_loader, 1):
            if device == "gpu":
                train_input = input.cuda()
                train_true = mask.cuda()
            else:
                train_input = input
                train_true = mask

            optimizer.zero_grad()

            train_pred = net(train_input)
            # output = torch.sigmoid(output)
            # output = output.cpu().detach()

            loss = criterion(train_pred, train_true)

            loss.backward()
            optimizer.step()

            iters += 1

            train_losses.append(loss.item())

            if output_dice:
                # # STEP:calculate dice for each class
                classes_dice = multi_class_dice(train_pred, train_true, dataset.num_classes, device="cpu")

                for i in range(1, dataset.num_classes):
                    cur_class_dice = classes_dice[i - 1]
                    epoch_classes_dice[i].append(cur_class_dice)

                # # STEP: calculate mean dice in this round
                mean_dice = sum(classes_dice) / len(classes_dice)

                # STEP: add into global
                epoch_mean_dice += mean_dice

        if output_dice:
            epoch_mean_dice = epoch_mean_dice / (len(train_loader) / 1)
            epoch_mean_classes_dice = [sum(epoch_classes_dice[i]) / len(epoch_classes_dice[i]) for i in range(1, dataset.num_classes)]

            train_loss = np.average(train_losses)

            et = datetime.datetime.now()

            log_file.write('epoch\t{}\ttrain_loss\t{:.4}\ttrain_mean_dice\t{:.4}\tREF\t{:.4}\tINS\t{:.4}\tDEL\t{:.4}\tINV\t{:.4}\tDUP\t{:.4}\tinvDUP\t{:.4}\ttime\t{}\n'.format(epoch, train_loss, epoch_mean_dice, epoch_mean_classes_dice[0], epoch_mean_classes_dice[1], epoch_mean_classes_dice[2], epoch_mean_classes_dice[3], epoch_mean_classes_dice[4], epoch_mean_classes_dice[5], (et - st).seconds))
            print('epoch {} - train_loss: {:.4} - train_mean_dice: {:.4} - REF: {:.4} - INS: {:.4} - DEL: {:.4} - INV: {:.4} - DUP: {:.4} - invDUP: {:.4} - time: {}'.format(epoch, train_loss, epoch_mean_dice, epoch_mean_classes_dice[0], epoch_mean_classes_dice[1], epoch_mean_classes_dice[2], epoch_mean_classes_dice[3], epoch_mean_classes_dice[4], epoch_mean_classes_dice[5], (et - st).seconds))
        else:
            train_loss = np.average(train_losses)

            et = datetime.datetime.now()
            log_file.write('epoch\t{}\ttrain_loss\t{:.4}\ttime\t{}\n'.format(epoch, train_loss, (et - st).seconds))
            print('epoch {} - train_loss: {:.4} - time: {}'.format(epoch, train_loss, (et - st).seconds))

        # # STEP: run validation
        val_classes_dice = [[] for i in range(dataset.num_classes)]
        val_mean_dice = 0

        st = datetime.datetime.now()
        net.eval().cuda()
        for val_batch, ((input, mask), file_name) in enumerate(val_loader, 1):
            if device == "gpu":
                val_input = input.cuda()
                val_true = mask.cuda()
            else:
                val_input = input
                val_true = mask

            val_pred = net(val_input)
            # pred = torch.sigmoid(pred)
            # pred = pred.cpu().detach()

            val_loss = criterion(val_pred, val_true)

            val_losses.append(val_loss.item())

            if output_dice:
                classes_dice = multi_class_dice(val_pred, val_true, dataset.num_classes, device="cpu")

                for i in range(1, dataset.num_classes):
                    cur_class_dice = classes_dice[i - 1]
                    val_classes_dice[i].append(cur_class_dice)

                # # STEP: calculate mean dice in this round
                mean_dice = sum(classes_dice) / len(classes_dice)

                val_mean_dice += mean_dice
        if output_dice:
            val_loss = np.average(val_losses)

            val_mean_dice = val_mean_dice / (len(val_loader) / 1)
            val_mean_classes_dice = [sum(val_classes_dice[i]) / len(val_classes_dice[i]) for i in range(1, dataset.num_classes)]

            et = datetime.datetime.now()
            print('epoch {} - val_loss: {:.4} - val_mean_dice: {:.4} - REF: {:.4} - INS: {:.4} - DEL: {:.4} - INV: {:.4} - DUP: {:.4} - invDUP: {:.4} - time: {}'.format(epoch, val_loss, val_mean_dice, val_mean_classes_dice[0], val_mean_classes_dice[1], val_mean_classes_dice[2], val_mean_classes_dice[3], val_mean_classes_dice[4], val_mean_classes_dice[5], (et - st).seconds))
            log_file.write('epoch\t{}\tval_loss\t{:.4}\tval_mean_dice\t{:.4}\tREF\t{:.4}\tINS\t{:.4}\tDEL\t{:.4}\tINV\t{:.4}\tDUP\t{:.4}\tinvDUP\t{:.4}\ttime\t{}\n'.format(epoch, val_loss, val_mean_dice, val_mean_classes_dice[0], val_mean_classes_dice[1], val_mean_classes_dice[2], val_mean_classes_dice[3], val_mean_classes_dice[4], val_mean_classes_dice[5], (et - st).seconds))
        else:
            val_loss = np.average(val_losses)

            et = datetime.datetime.now()

            print('epoch {} - val_loss: {:.4} - time: {}'.format(epoch, val_loss, (et - st).seconds))
            log_file.write('epoch\t{}\tval_loss\t{:.4}\ttime\t{}\n'.format(epoch, val_loss,  (et - st).seconds))

        #print(optimizer.param_groups[0]['lr'])

        #if val_mean_dice >= 0.97:
        #    optimizer.param_groups[0]['lr'] = 1e-5

        #if val_mean_dice >= 0.985:
        #    optimizer.param_groups[0]['lr'] = 1e-6
        # print(optimizer.param_groups[0]['lr'])

        if scheduler is not None:
            scheduler.step()

        early_stopping(val_mean_dice, net, epoch)
        # print(early_stopping.early_stop, os.path.join(root_path, '{}.pth'.format(model_name)))
        # if early_stopping.early_stop or optimizer.param_groups[0]['lr'] < threshold_lr:

        if early_stopping.early_stop:
            print("Early stopping")
            # 结束模型训练
            break

    print('----------------------------------------------------------')
    print('save epoch {}'.format(early_stopping.save_epoch))
    print('stoped epoch {}'.format(epoch))
    print('----------------------------------------------------------')


def main():
    # # STEP: define data transform
    # image_norm_transform = image_transforms.ImageNormalize()
    image_norm_transform = None

    center_crop_transform = image_transforms.CenterCrop(crop_size)
    img2tensor_transform = image_transforms.NpyToTensor()
    mask2tensor_transform = image_transforms.MaskToTensor()

    # # STEP: load model and val set
    train_set = dataset.Dataset(train_path, 'model', fold, image_norm_transform=image_norm_transform, img2tensor_transform=img2tensor_transform, center_crop_transform=center_crop_transform, mask2tensor_transform=mask2tensor_transform)
    train_loader = DataLoader(train_set, batch_size=batch_size, shuffle=True, num_workers=12, drop_last=True)

    val_set = dataset.Dataset(val_path, 'val', fold, image_norm_transform=image_norm_transform, img2tensor_transform=img2tensor_transform, center_crop_transform=center_crop_transform, mask2tensor_transform=mask2tensor_transform)
    val_loader = DataLoader(val_set, batch_size=batch_size, shuffle=False, num_workers=12, drop_last=True)

    # # STEP: define the network
    if model == "unet":
        depth = 2  # unet conv depth
        loss_name = 'ce'  # dice, bce, wbce, dual, wdual
        optimizer_type = 'adam'  # adam, sgd, RMSprop
        scheduler_type = 'no'  # ReduceLR, StepLR, poly

        if device == "gpu":
            net = Unet(num_classes=dataset.num_classes, depth=depth).cuda()
        else:
            net = Unet(num_classes=dataset.num_classes, depth=depth)
    elif model == "liteunet":
        depth = 2  # unet conv depth
        loss_name = 'ce'  # dice, bce, wbce, dual, wdual
        optimizer_type = 'adam'  # adam, sgd, RMSprop
        scheduler_type = 'no'  # ReduceLR, StepLR, poly

        if device == "gpu":
            net = LiteUnet(num_classes=dataset.num_classes, depth=depth).cuda()
        else:
            net = LiteUnet(num_classes=dataset.num_classes, depth=depth)

    elif model == "miniunet":
        depth = 2  # unet conv depth
        loss_name = 'ce'  # dice, bce, wbce, dual, wdual
        optimizer_type = 'adam'  # adam, sgd, RMSprop
        scheduler_type = 'no'  # ReduceLR, StepLR, poly

        if device == "gpu":
            net = MiniUnet(num_classes=dataset.num_classes, depth=depth).cuda()
        else:
            net = MiniUnet(num_classes=dataset.num_classes, depth=depth)

    elif model == "fcn":
        loss_name = 'ce'  # dice, bce, wbce, dual, wdual
        optimizer_type = 'rmsprop'  # adam, sgd, rmsprop
        scheduler_type = 'no'  # ReduceLR, StepLR, poly

        if device == "gpu":
            vgg_model = VGGNet(requires_grad=True, remove_fc=True).cuda()
            net = FCNs(pretrained_net=vgg_model, n_class=dataset.num_classes).cuda()
        else:
            vgg_model = VGGNet(requires_grad=True, remove_fc=True)
            net = FCNs(pretrained_net=vgg_model, n_class=dataset.num_classes)

    elif model == "deeplabv3":
        loss_name = 'ce'  # dice, bce, wbce, dual, wdual
        optimizer_type = 'rmsprop'  # adam, sgd, RMSprop
        scheduler_type = 'no'  # ReduceLR, StepLR, poly

        if device == "gpu":
            net = DeepLabV3(num_classes=dataset.num_classes).cuda()
        else:
            net = DeepLabV3(num_classes=dataset.num_classes)

    else:
        print("[No this model]")
        exit()

    # # STEP: define loss function
    if loss_name == 'dice':
        if device == "gpu":
            criterion = Dice_Loss(dataset.num_classes).cuda()
        else:
            criterion = Dice_Loss(dataset.num_classes)

    elif loss_name == "ce":
        if device == "gpu":
            criterion = CE_Dice_Loss(dataset.num_classes).cuda()
        else:
            criterion = CE_Dice_Loss(dataset.num_classes)
    else:
        criterion = None
        print("No this criterion type: {}".format(loss_name))
        exit()

    # # STEP: define early stop
    early_stopping = EarlyStopping(early_stop_patience, verbose=True, delta=early_stop__eps, path=os.path.join(train_path, '{}.pth'.format(model_name)))

    # # STEP: define optimizer
    if optimizer_type == 'adam':
        optimizer = torch.optim.Adam(net.parameters(), lr=initial_lr, weight_decay=weight_decay)
    elif optimizer_type == "rmsprop":
        optimizer = torch.optim.RMSprop(net.parameters(), lr=initial_lr, momentum=0.999, weight_decay=1e-8)
    elif optimizer_type == "sgd":
        optimizer = torch.optim.SGD(net.parameters(), lr=0.1, momentum=0.9)
    else:
        optimizer = None
        print("No this optimizer type: {}".format(optimizer_type))
        exit()

    # # STEP: define learning rate
    if scheduler_type == 'StepLR':
        scheduler = lr_scheduler.StepLR(optimizer, step_size=40, gamma=0.5)
    elif scheduler_type == 'ReduceLR':
        scheduler = lr_scheduler.ReduceLROnPlateau(optimizer, mode='max', patience=5)
    elif scheduler_type == 'poly':
        scheduler = PolyLR(optimizer, max_iter=n_epoch, power=0.9)
    else:
        scheduler = None

    # # STEP: begin training
    train(train_loader, val_loader, net, criterion, optimizer, scheduler, early_stopping, n_epoch, 0)


if __name__ == '__main__':
    main()

    log_file.close()
