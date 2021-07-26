Print *, "Initialising grid"
            call conf%make_mesh(ierr)
            allocate (x_pass_through(size(conf%xx), 2))
            x_pass_through(:, 1) = reshape(conf%xx, (/size(conf%xx)/))
            x_pass_through(:, 2) = reshape(conf%yy, (/size(conf%yy)/))
            Laplacian_COO = laplacian2d(conf%x, conf%y, conf%BCx, conf%BCy, ierr)
            Laplacian_CSR = coo_to_csr(Laplacian_COO, ierr)