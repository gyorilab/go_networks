from go_networks.generate import get_genes_for_go_id_direct, get_genes_for_go_id


def test_get_genes_for_go_id():
    go_id = 'GO:0061564'
    genes_direct = get_genes_for_go_id_direct(go_id)
    genes = get_genes_for_go_id(go_id)
    assert len(genes) > len(genes_direct)
    assert len(genes) == 523
