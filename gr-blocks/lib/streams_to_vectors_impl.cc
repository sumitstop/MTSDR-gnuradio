/* -*- c++ -*- */
/*
 * Copyright 2012 Free Software Foundation, Inc.
 *
 * This file is part of GNU Radio
 *
 * GNU Radio is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 *
 * GNU Radio is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GNU Radio; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street,
 * Boston, MA 02110-1301, USA.
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "streams_to_vectors_impl.h"
#include <gnuradio/io_signature.h>

namespace gr {
  namespace blocks {

    streams_to_vectors::sptr streams_to_vectors::make(size_t itemsize, size_t nitems_per_block)
    {
      return gnuradio::get_initial_sptr(new streams_to_vectors_impl(itemsize, nitems_per_block));
    }

    streams_to_vectors_impl::streams_to_vectors_impl(size_t itemsize, size_t nitems_per_block)
      : sync_decimator ("streams_to_vectors",
			   io_signature::make2 (2, 2, itemsize, itemsize),
			   io_signature::make2 (2, 2, itemsize * nitems_per_block, itemsize * nitems_per_block),
			   nitems_per_block)
    {
    }

    int
    streams_to_vectors_impl::work(int noutput_items,
				 gr_vector_const_void_star &input_items,
				 gr_vector_void_star &output_items)
    {
      size_t block_size = output_signature()->sizeof_stream_item (0);
      size_t block_size_1 = output_signature()->sizeof_stream_item (0);

      const char *in = (const char *) input_items[0];
      const char *in_1 = (const char *) input_items[1];
      char *out = (char *) output_items[0];
      char *out_1 = (char *) output_items[1];

      memcpy (out, in, noutput_items * block_size);
      memcpy (out_1, in_1, noutput_items * block_size_1);

      return noutput_items;
    }

  } /* namespace blocks */
} /* namespace gr */
