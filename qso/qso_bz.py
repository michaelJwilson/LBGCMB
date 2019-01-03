import  numpy  as  np


def qso_bz(z):
    return 0.278 * ((1. + z)**2. - 6.565) + 2.393


if __name__ == "__main__":
    print("\n\nWelcome to QSO b(z).\n\n")

    print  qso_bz(1.16)
    print  qso_bz(2.61)

    print("\n\nDone.\n\n")
