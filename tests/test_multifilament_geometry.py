from apexloss.multifilament_fem import hex19_centers, secondary304_centers


def test_topologies():
    assert len(hex19_centers()) == 19
    c, g = secondary304_centers()
    assert len(c) == 304
    assert len(set(g.tolist())) == 16
