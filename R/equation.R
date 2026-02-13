#' @importFrom generics equation
#' @export
generics::equation

#' Extract model equations
#'
#' @param object An `affectOU` object
#' @param type Output format; one of `"plain"`, `"expression"`, `"latex"`, or `"code"`
#' @param inline If TRUE, insert numeric parameter values into equations.
#'   If FALSE, keep symbolic notation and define parameters separately.
#' @param digits Number of digits for rounding numeric values
#' @param ... Additional arguments (unused)
#'
#' @export
#'
#' @examples
#' # Plain text equation for 1D model
#' model <- affectOU()
#' cat(equation(model))
#'
#' # Inlined Latex equation for 2D model
#' model <- affectOU(ndim = 2)
#' cat(equation(model, type = "latex", inline = TRUE))
#'
#' # Expression output for 1D model
#' model <- affectOU()
#' equation(model, type = "expression")
#'
#' # R code output for 1D model
#' model <- affectOU()
#' cat(equation(model, type = "code"))
#'
equation.affectOU <- function(object, type = c("plain", "latex", "expression", "code"),
                              inline = FALSE, digits = 3, ...) {
  type <- match.arg(type)

  if (object[["ndim"]] == 1) {
    # Ensure parameters are numeric scalars
    for (param in c("theta", "mu", "gamma", "sigma")) {
      object[["parameters"]][[param]] <- as.numeric(object[["parameters"]][[param]])
    }

    equation_ou_1d(object, type, inline, digits)
  } else {
    equation_ou_nd(object, type, inline, digits)
  }
}


# 1D equations ---------------------------------------------------------------

equation_ou_1d <- function(object, type, inline, digits) {
  theta <- round(object$parameters$theta, digits)
  mu <- round(object$parameters$mu, digits)
  gamma <- round(object$parameters$gamma, digits)
  sigma <- round(object$parameters$sigma, digits)
  switch(type,
    plain = equation_ou_1d_plain(theta, mu, gamma, sigma, inline),
    latex = equation_ou_1d_latex(theta, mu, gamma, sigma, inline),
    expression = equation_ou_1d_expression(theta, mu, gamma, sigma, inline),
    code = equation_ou_1d_code(theta, mu, gamma, sigma, inline)
  )
}

equation_ou_1d_plain <- function(theta, mu, gamma, sigma, inline) {
  if (inline) {
    sprintf("dX(t) = %s * (%s - X(t)) dt + %s dW(t)", theta, mu, gamma)
  } else {
    paste0(
      "dX(t) = theta * (mu - X(t)) dt + gamma dW(t)\n\n",
      "where:\n",
      sprintf("  theta = %s\n", theta),
      sprintf("  mu    = %s\n", mu),
      sprintf("  gamma = %s\n", gamma),
      sprintf("  sigma = %s\n", sigma)
    )
  }
}

equation_ou_1d_latex <- function(theta, mu, gamma, sigma, inline) {
  if (inline) {
    sprintf(
      "dX(t) = %s \\left( %s - X(t) \\right) dt + %s \\, dW(t)",
      theta, mu, gamma
    )
  } else {
    paste0(
      "dX(t) = \\theta \\left( \\mu - X(t) \\right) dt + \\gamma \\, dW(t)\n\n",
      "\\text{where:}\n",
      "\\begin{align*}\n",
      sprintf("  \\theta &= %s \\\\\n", theta),
      sprintf("  \\mu &= %s \\\\\n", mu),
      sprintf("  \\gamma &= %s \\\\\n", gamma),
      sprintf("  \\sigma &= %s\n", sigma),
      "\\end{align*}\n"
    )
  }
}

equation_ou_1d_expression <- function(theta, mu, gamma, sigma, inline) {
  if (inline) {
    bquote(dX(t) == .(theta) * (.(mu) - X(t)) * dt + .(gamma) * dW(t))
  } else {
    list(
      equation = quote(dX(t) == theta * (mu - X(t)) * dt + gamma * dW(t)),
      theta = theta,
      mu = mu,
      gamma = gamma,
      sigma = sigma
    )
  }
}

equation_ou_1d_code <- function(theta, mu, gamma, sigma, inline) {
  if (inline) {
    sprintf("dX <- %s * (%s - X) * dt + %s * dW", theta, mu, gamma)
  } else {
    paste0(
      sprintf("theta <- %s\n", theta),
      sprintf("mu <- %s\n", mu),
      sprintf("gamma <- %s\n", gamma),
      sprintf("sigma <- %s\n", sigma),
      "\n",
      "dX <- theta * (mu - X) * dt + gamma * dW\n"
    )
  }
}


# Multivariate equations -----------------------------------------------------

equation_ou_nd <- function(object, type, inline, digits) {
  theta <- round(object$parameters$theta, digits)
  mu <- round(object$parameters$mu, digits)
  gamma <- round(object$parameters$gamma, digits)
  sigma <- round(gamma %*% t(gamma), digits)

  switch(type,
    plain = equation_ou_nd_plain(theta, mu, gamma, sigma, inline),
    latex = equation_ou_nd_latex(theta, mu, gamma, sigma, inline),
    expression = equation_ou_nd_expression(theta, mu, gamma, sigma, inline),
    code = equation_ou_nd_code(theta, mu, gamma, sigma, inline)
  )
}

equation_ou_nd_plain <- function(theta, mu, gamma, sigma, inline) {
  if (inline) {
    paste0(
      "dX(t) = Theta (mu - X(t)) dt + Gamma dW(t)\n\n",
      "Theta =\n", format_matrix_plain(theta), "\n\n",
      "mu = [", paste(mu, collapse = ", "), "]\n\n",
      "Gamma =\n", format_matrix_plain(gamma), "\n\n",
      "Sigma =\n", format_matrix_plain(sigma), "\n"
    )
  } else {
    paste0(
      "dX(t) = Theta (mu - X(t)) dt + Gamma dW(t)\n\n",
      "where:\n\n",
      "Theta =\n", format_matrix_plain(theta), "\n\n",
      "mu = [", paste(mu, collapse = ", "), "]\n\n",
      "Gamma =\n", format_matrix_plain(gamma), "\n\n",
      "Sigma =\n", format_matrix_plain(sigma), "\n"
    )
  }
}

equation_ou_nd_latex <- function(theta, mu, gamma, sigma, inline) {
  if (inline) {
    paste0(
      "d\\mathbf{X}(t) = ",
      format_matrix_latex(theta),
      " \\left( ",
      format_vector_latex(mu),
      " - \\mathbf{X}(t) \\right) dt + ",
      format_matrix_latex(gamma),
      " \\, d\\mathbf{W}(t)\n"
    )
  } else {
    paste0(
      "d\\mathbf{X}(t) = \\boldsymbol{\\Theta} ",
      "\\left( \\boldsymbol{\\mu} - \\mathbf{X}(t) \\right) dt + ",
      "\\boldsymbol{\\Gamma} \\, d\\mathbf{W}(t)\n\n",
      "\\text{where:}\n",
      "\\begin{align*}\n",
      "  \\boldsymbol{\\Theta} &= ", format_matrix_latex(theta), " \\\\\n",
      "  \\boldsymbol{\\mu} &= ", format_vector_latex(mu), " \\\\\n",
      "  \\boldsymbol{\\Gamma} &= ", format_matrix_latex(gamma), " \\\\\n",
      "  \\boldsymbol{\\Sigma} &= ", format_matrix_latex(sigma), "\n",
      "\\end{align*}\n"
    )
  }
}

equation_ou_nd_expression <- function(theta, mu, gamma, sigma, inline) {
  # Expressions with matrices are impractical, so always return a list
  list(
    equation = quote(d * bold(X)(t) == bold(Theta) * (bold(mu) - bold(X)(t)) * dt + bold(Gamma) * d * bold(W)(t)),
    theta = theta,
    mu = mu,
    gamma = gamma,
    sigma = sigma
  )
}

equation_ou_nd_code <- function(theta, mu, gamma, sigma, inline) {
  ndim <- nrow(theta)
  paste0(
    sprintf("Theta <- %s\n", format_matrix_code(theta)),
    sprintf("mu <- %s\n", format_vector_code(mu)),
    sprintf("Gamma <- %s\n", format_matrix_code(gamma)),
    sprintf("Sigma <- %s\n", format_matrix_code(sigma)),
    "\n",
    "dX <- Theta %*% (mu - X) * dt + Gamma %*% dW\n"
  )
}


# Formatting helpers ---------------------------------------------------------

format_matrix_plain <- function(m) {
  # Use capture.output with print to get clean matrix representation
  lines <- utils::capture.output(print(m, quote = FALSE))
  paste(lines, collapse = "\n")
}

format_matrix_latex <- function(m) {
  rows <- apply(m, 1, function(row) paste(row, collapse = " & "))
  paste0(
    "\\begin{pmatrix}\n  ",
    paste(rows, collapse = " \\\\\n  "),
    "\n\\end{pmatrix}"
  )
}

format_vector_latex <- function(v) {
  paste0(
    "\\begin{pmatrix} ",
    paste(v, collapse = " \\\\ "),
    " \\end{pmatrix}"
  )
}

format_matrix_code <- function(m) {
  values <- paste(as.vector(m), collapse = ", ")
  sprintf("matrix(c(%s), nrow = %d, ncol = %d)", values, nrow(m), ncol(m))
}

format_vector_code <- function(v) {
  sprintf("c(%s)", paste(v, collapse = ", "))
}
