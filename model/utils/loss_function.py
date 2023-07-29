
from torch.nn.modules.loss import _Loss
import torch.nn.functional as F
import torch.nn as nn


def dice_coeff(y_pred, y_true, smooth=1):
    assert y_pred.size() == y_true.size()

    y_pred = y_pred[:, 0].contiguous().view(-1)

    y_true = y_true[:, 0].contiguous().view(-1)

    # print(y_pred.device, y_true.device)
    intersection = (y_pred * y_true).sum()

    dsc = (2. * intersection + smooth) / (
            y_pred.sum() + y_true.sum() + smooth
    )

    # print(y_pred, y_pred.sum(), np.shape(y_pred))
    # print(y_true, y_true.sum(), np.shape(y_true))
    # print(00000, np.shape(np.where(y_pred > 0.1)))
    # print(00000, np.shape(np.where(y_true > 0.1)))
    #
    # print(intersection)
    # print(dsc)

    return dsc


def multi_class_dice(y_pred, y_true, num_classes, device='gpu'):
    """
    calculate dice value for reach class
    """

    y_pred = F.softmax(y_pred, dim=1).float()
    # soft_max = nn.Softmax(dim=1).cuda()
    # y_pred = soft_max(y_pred)

    class_dice = []

    # start with 1 to except backgrpund
    for i in range(1, num_classes):
        if device == "cpu":
            class_dice.append(dice_coeff(y_pred[:, i:i + 1, :], y_true[:, i:i + 1, :]).cpu().item())
        else:
            class_dice.append(dice_coeff(y_pred[:, i:i + 1, :], y_true[:, i:i + 1, :]))

    return class_dice


class Dice_Loss(_Loss):

    def __init__(self, num_classes):
        super(Dice_Loss, self).__init__()
        self.num_classes = num_classes

    def forward(self, y_pred, y_true):

        class_dice = multi_class_dice(y_pred, y_true, self.num_classes)

        mean_dice = sum(class_dice) / len(class_dice)

        return 1 - mean_dice


# class BCE_Dice_Loss(_Loss):
#     def __init__(self, num_classes, smooth=0, weight=[1.0, 1.0]):
#         super(BCE_Dice_Loss, self).__init__()
#         self.bce_loss = nn.BCEWithLogitsLoss().cuda()
#         self.dice_loss = Dice_Loss(num_classes=num_classes).cuda()
#         self.weight = weight
#         self.smooth = smooth
#         self.num_classes = num_classes
#
#     def forward(self, inputs, targets):
#         return self.weight[0] * self.bce_loss(inputs, targets * (1 - self.smooth) + self.smooth / self.num_classes) + self.weight[1] * self.dice_loss(inputs, targets)
#

class CE_Dice_Loss(_Loss):
    def __init__(self, num_classes, smooth=0, weight=[1.0, 1.0]):
        super(CE_Dice_Loss, self).__init__()
        self.bce_loss = nn.CrossEntropyLoss().cuda()
        self.dice_loss = Dice_Loss(num_classes=num_classes).cuda()
        self.weight = weight
        self.smooth = smooth
        self.num_classes = num_classes

    def forward(self, inputs, targets):
        return self.weight[0] * self.bce_loss(inputs, targets * (1 - self.smooth) + self.smooth / self.num_classes) + self.weight[1] * self.dice_loss(inputs, targets)


