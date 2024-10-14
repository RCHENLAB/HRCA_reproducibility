import torch

from utils import common_functions as c_f
from losses.mixins import WeightRegularizerMixin
from losses.nca_loss import NCALoss


class ProxyNCALoss(WeightRegularizerMixin, NCALoss):
    def __init__(self, num_classes, embedding_size, **kwargs):
        super().__init__(**kwargs)
        self.proxies = torch.nn.Parameter(torch.Tensor(num_classes, embedding_size))
        self.weight_init_func(self.proxies)
        self.proxy_labels = torch.arange(num_classes)
        self.add_to_recordable_attributes(list_of_names=["num_classes"], is_stat=False)

    def cast_types(self, dtype, device):
        self.proxies.data = c_f.to_device(self.proxies.data, device=device, dtype=dtype)

    def compute_loss(self, embeddings, labels, indices_tuple):
        dtype, device = embeddings.dtype, embeddings.device
        self.cast_types(dtype, device)
        loss_dict = self.nca_computation(
            embeddings,
            self.proxies,
            labels,
            c_f.to_device(self.proxy_labels, labels),
            indices_tuple,
        )
        self.add_weight_regularization_to_loss_dict(loss_dict, self.proxies)
        return loss_dict
