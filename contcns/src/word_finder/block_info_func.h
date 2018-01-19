/// block info function
typedef int (*BlkInfoFunc)(idx soff, idx block_size);

int
general_bid_func(idx soff, idx block_size)
{
    return (int)(soff / block_size);
}

int
general_boff_func(idx soff, idx block_size)
{
    return (int)(soff % block_size);
}

int
fast_bid_func(idx soff, idx)
{
    return (int)(soff >> sBidShift);
}

int 
fast_boff_func(idx soff, idx)
{
    return (int)(soff & sBoffMask);
}

typedef idx (*SoffFunc)(idx bid, idx boff, idx block_size);

idx
general_soff_func(idx bid, idx boff, idx block_size)
{
	idx r = bid * block_size + boff;
	return r;
}

idx
fast_soff_func(idx bid, idx boff, idx block_size)
{
	idx r = (bid << sBidShift) | boff;
	return r;
}

void
decide_block_info_func(const idx n)
{
    if (n & (n-1)) { /// not power of 2
        BidFunc = general_bid_func;
        BoffFunc = general_boff_func;
		sSoffFunc = general_soff_func;
    } else {
        idx m = 1;
        while (m != n) {
            ++sBidShift;
            m *= 2;
        }
        sBoffMask = (U64_ONE << sBidShift) - 1;

        BidFunc = fast_bid_func;
        BoffFunc = fast_boff_func;
		sSoffFunc = fast_soff_func;
    }
}

