a =\
[
    [],
    [1/2.],
    [0,     1/2.],
    [0,     0,      1.]
]
b = [1/6., 1/3., 1/3., 1/6.]
c = [0, 1/2., 1/2., 1]
steps = len(b)


def appendto(_list, _element):
    _list.append(_element)
    return _list


def rungekutta(t, y, f, h):
    k = []
    for i in xrange(steps):
        k.append(h * f(t + c[i]*h, y + reduce(lambda x, y: x+y, map(lambda x, y: x*y, a[i], k), 0)))
    return y + reduce(lambda x, y: x+y, map(lambda x, y: x*y, b, k))


def rungekuttavec(t, yv, fv, h):
    kv = []
    for i in xrange(steps):
        # print "step = " + str(i)
        # print kv
        # print a[i]
        ki = []
        for idx in xrange(len(yv)):
            # print "idx = " + str(idx)
            tyv = map(lambda _factor, _vector: map(lambda _vector_elem: _factor * _vector_elem, _vector), a[i], kv)
            tyv = reduce(lambda _vec1, _vec2: map(
                lambda _vec1_elem, _vec2_elem: _vec1_elem + _vec2_elem,
                _vec1, _vec2
            ), tyv, [0. for q in xrange(len(yv))])
            tyv = map(lambda _vec1_elem, _vec2_elem: _vec1_elem + _vec2_elem, tyv, yv)
            ki.append(h * fv[idx](t + c[i]*h, tyv))
        kv.append(ki)

    res = map(lambda _factor, _vector: map(lambda _vector_elem: _factor * _vector_elem, _vector), b, kv)
    res = reduce(lambda _vec1, _vec2: map(
        lambda _vec1_elem, _vec2_elem: _vec1_elem + _vec2_elem,
        _vec1, _vec2
    ), res, [0. for q in xrange(len(yv))])
    res = map(lambda _vec1, _vec2: _vec1 + _vec2, yv, res)
    return res


def solve(fv, y0v, step, stepcount = 100, x0 = 0.0):
    x = x0
    yv = y0v
    stepcount = int(stepcount)
    graph = [[] for i in xrange(len(fv))]
    for i in xrange(stepcount):
        yv = rungekuttavec(x, yv, fv, step)
        x += step
        for j in xrange(len(fv)):
            graph[j].append([x, yv[j]])
    return graph