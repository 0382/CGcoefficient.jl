function test_xGaunt(test_range::AbstractArray)
    for l1 in test_range, l2 in test_range, l3 in test_range
        for m1 in -l1:l1, m2 in -l2:l2, m3 in -l3:l3
            g = xGaunt(l1, l2, l3, m1, m2, m3)
            w = exact_sqrt((2l1 + 1) * (2l2 + 1) * (2l3 + 1)//4) * threeJ(l1, l2, l3, m1, m2, m3) * threeJ(l1, l2, l3, 0, 0, 0)
            @test g == w
        end
    end
end